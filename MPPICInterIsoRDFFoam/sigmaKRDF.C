/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "sigmaKRDF.H"
#include "zoneDistribute.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::sigmaKRDF::sigmaKRDF
(
    const fvMesh& mesh,
    const immiscibleIncompressibleTwoPhaseMixture& mixture,
    //const isoAdvection& advector,
    isoAdvection& advector,
    const dictionary& dict
)
:
    mesh_(mesh),
    mixture_(mixture),
    advector_(advector),
    rdfPtr_(new reconstructedDistanceFunction(mesh)),
    bandWidth_(dict.getOrDefault<label>("RDFbandWidth", 4)),
    nSmoothIter_(dict.getOrDefault<label>("RDFnSmoothIter", 0)),
    deltaN_
    (
        "deltaN",
        dimless,
        // v2 清理点 #4：数值跟网格缩放
        //   特征长度 ≈ cbrt(平均 cell 体积)；1e-8 除特征长度得到 dimless 稳定项
        1e-8 / Foam::cbrt(Foam::gAverage(mesh.V()))
    ),
    sigma_
    (
        "sigma",
        // N/m = kg/s²
        dimensionSet(1, 0, -2, 0, 0, 0, 0),
        dict.get<scalar>("RDFsigma")
    )
{
    Info<< "sigmaKRDF initialized:" << nl
        << "    bandWidth    = " << bandWidth_ << nl
        << "    nSmoothIter  = " << nSmoothIter_ << nl
        << "    deltaN       = " << deltaN_.value() << " (dimless)" << nl
        << "    sigma        = " << sigma_.value() << " N/m" << endl;
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

void Foam::sigmaKRDF::update()
{
    // === share 关键点：不调 surf().reconstruct()，复用主 advector 的 PLIC 结果 ===
    // advect() 在 alphaEqnSubCycle.H 内部已经重构过，此时 surf() 就是最新 α 对应的 PLIC 平面。

    // 标记距界面 ≤ bandWidth 层的 cell（窄带）
    rdfPtr_->markCellsNearSurf
    (
        advector_.surf().interfaceCell(),
        bandWidth_
    );

    // 并行窄带通信（跨 rank 边界传递 PLIC 平面信息）
    zoneDistribute& exchange = zoneDistribute::New(mesh_);

    // 构造 RDF：从每个窄带 cell 找最近的 PLIC 平面，RDF = 该 cell 中心到该平面的带符号距离
    rdfPtr_->constructRDF
    (
        rdfPtr_->nextToInterface(),
        advector_.surf().centre(),
        advector_.surf().normal(),
        exchange,
        true    // updateStencil
    );

    rdfPtr_->correctBoundaryConditions();
}


Foam::tmp<Foam::volScalarField> Foam::sigmaKRDF::sigmaK() const
{
    const reconstructedDistanceFunction& rdf = rdfPtr_();

    // --- κ_RDF = -div(grad(RDF)/|grad(RDF)|) ---------------------------------
    // 量纲：grad(RDF) = grad(dimLength) = dimless（这是之前踩过的坑点，保留）
    volVectorField gRDF(fvc::grad(rdf));
    volVectorField nRDF(gRDF / (mag(gRDF) + deltaN_));
    nRDF.rename("nRDF");                    // ★ 新增：强制字段名为 "nRDF"，让 divScheme 查 div(nRDF)
    nRDF.correctBoundaryConditions();

    // 可选法向平滑（v2 清理点 #2：默认 0，不跑；
    // 若未来开启，是平滑 n̂ 而非 κ，不会抹平三相接触线尖锐曲率）
    for (label i = 0; i < nSmoothIter_; ++i)
    {
        nRDF = fvc::average(fvc::interpolate(nRDF));
        nRDF.correctBoundaryConditions();
    }

    volScalarField kRDF(-fvc::div(nRDF));

    // --- 连续 mask（v2 清理点 #3：用 cellDistLevel 线性渐变） ----------------
    // cellDistLevel ∈ {0, 1, 2, ..., bandWidth} 表示距界面的 cell 层数
    //   level = 0（界面 cell）        → mask = 1  → 纯 RDF
    //   level = bandWidth（窄带边界） → mask = 0  → 纯 CSF
    //   level 超过 bandWidth（窄带外）→ min 会给出 ≤0，max 截断到 0 → 纯 CSF
    volScalarField mask
    (
        max
        (
            scalar(0),
            min
            (
                scalar(1),
                (scalar(bandWidth_) - rdf.cellDistLevel()) / scalar(bandWidth_)
            )
        )
    );

    // --- κ 层 blend，外乘 σ（v2 清理点 #5） ----------------------------------
    // 数学等价于旧版的 mask*(σ*κ_RDF) + (1-mask)*mixture.sigmaK()，
    // 但不需要手工解析 σ，更 OpenFOAM-native。
    return mask*sigma_*kRDF + (scalar(1) - mask)*mixture_.sigmaK();
}


// ************************************************************************* //
