// ============================================================================
// alphaSuSp.H （方案 A 本地自包含版）
// ----------------------------------------------------------------------------
// 功能：初始化 α 方程的 Internal 源项 Sp、Su，并吸收 fvOptions 对 alpha1 的
// 线性/显式贡献。isoAdvector::advect(Sp, Su) 接受这两个 Internal 字段。
//
// 与标准 interFoam/alphaSuSp.H 的差别：
//   - 不计算 divU / phiCN：isoAdvector 几何输运内部自洽，外部不需要提供这些。
//   - 其余 fvOptions 吸收逻辑一致，保证 fvOptions 兼容性。
// ============================================================================

volScalarField::Internal Sp
(
    IOobject
    (
        "Sp",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dimless/dimTime, Zero)
);

volScalarField::Internal Su
(
    IOobject
    (
        "Su",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dimless/dimTime, Zero)
);

// 吸收 fvOptions 提供的 alpha1 源项
if (fvOptions.appliesToField(alpha1.name()))
{
    tmp<fvScalarMatrix> tSusp(fvOptions(alpha1));
    const fvScalarMatrix& Susp = tSusp();

    const scalarField& d = Susp.diag();
    const scalarField& src = Susp.source();
    const scalarField& V = mesh.V();

    forAll(Sp, i)
    {
        Sp[i] += d[i]/V[i];
        Su[i] -= src[i]/V[i];
    }
}
