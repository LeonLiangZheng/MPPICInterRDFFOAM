#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

// Need to store mu as incompressibleTwoPhaseMixture does not store it
volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    mixture.mu(),
    calculatedFvPatchScalarField::typeName
);


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    mesh.solutionDict().subDict("PIMPLE"),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// alphac must be constructed before the cloud
// so that the drag-models can find it
volScalarField alphac
(
    IOobject
    (
        "alphac",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, Zero),
    zeroGradientFvPatchScalarField::typeName
);
alphac.oldTime();

volScalarField alphacRho(alphac*rho);
alphacRho.oldTime();

Info<< "Constructing kinematicCloud " << endl;
basicKinematicCloud kinematicCloud
(
    "kinematicCloud",
    rho,
    U,
    mu,
    g
);

// Particle fraction upper limit
scalar alphacMin
(
    1.0
  - (
        kinematicCloud.particleProperties().subDict("constantProperties")
       .get<scalar>("alphaMax")
    )
);

// Update alphac from the particle locations
alphac = max(1.0 - kinematicCloud.theta(), alphacMin);
alphac.correctBoundaryConditions();

surfaceScalarField alphacf("alphacf", fvc::interpolate(alphac));

// Phase mass flux
surfaceScalarField alphaRhoPhic("alphaRhoPhic", alphacf*rhoPhi);

// Volumetric phase flux
surfaceScalarField alphaPhic("alphaPhic", alphacf*phi);

autoPtr
<
    PhaseCompressibleTurbulenceModel
    <
        immiscibleIncompressibleTwoPhaseMixture
    >
>turbulence
(
    PhaseCompressibleTurbulenceModel
    <
        immiscibleIncompressibleTwoPhaseMixture
    >::New
    (
        alphac,
        rho,
        U,
        alphaRhoPhic,
        rhoPhi,
        mixture
    )
);

#include "createMRF.H"

// ============================================================================
//                         【方案 A 新增】构造 isoAdvector
// ============================================================================
// 参数仅 (alpha1, phi, U)：isoAdvector 内部会从 phi 做 PLIC 重构。
// 注意：这里传的是 α1（而非 αc·α1），因为我们在 alphaEqn.H 里对 α 方程做了
// 方案 A 变形，isoAdvector 解的是变形后的 ∂α1/∂t + ∇·(α1·U) = Su'+Sp'·α1，
// 其源项承载了 αc 的影响。
//
// fvSolution 需要在 solvers 字典里加 "alpha.water.*" 子字典，含
// nAlphaSubCycles、reconstructionScheme (plicRDF 或 isoAlpha) 等参数。
// ============================================================================
isoAdvection advector(alpha1, phi, U);


// ============================================================================
//                        【方案 B 新增】RDF 曲率切换
// ============================================================================
// fvSolution 的 PIMPLE 子字典可含：
//   curvatureModel   RDF;    // 或 CSF；默认 CSF（= MPPICInterIsoFoam 行为）
//   RDFbandWidth     4;
//   RDFnSmoothIter   0;
//
// CSF 模式下不构造 sigmaKRDF，零额外开销；UEqn/pEqn 走 mixture.sigmaK() 分支。
// ============================================================================
const dictionary& pimpleDict = pimple.dict();

const word curvatureModel
(
    pimpleDict.getOrDefault<word>("curvatureModel", "CSF")
);

const bool useRDFCurvature = (curvatureModel == "RDF");

Info<< "Curvature model: " << curvatureModel
    << (useRDFCurvature ? " (RDF in band, CSF outside)" : "")
    << endl;

autoPtr<sigmaKRDF> rdfCurvPtr;
if (useRDFCurvature)
{
    rdfCurvPtr.reset(new sigmaKRDF(mesh, mixture, advector, pimpleDict));
}
