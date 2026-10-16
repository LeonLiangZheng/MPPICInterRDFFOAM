/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | www.openfoam.com
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Part of MPPICInterRDFFoam. Not part of upstream OpenFOAM.
-------------------------------------------------------------------------------
License
    GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "sigmaKRDF.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvcAverage.H"
#include "fvcSurfaceIntegrate.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigmaKRDF::sigmaKRDF
(
    immiscibleIncompressibleTwoPhaseMixture& mixture,
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U,
    const dictionary& dict
)
:
    mesh_(alpha1.mesh()),
    mixture_(mixture),
    alpha1_(alpha1),
    phi_(phi),
    U_(U),
    advectorPtr_(nullptr),
    rdfPtr_(nullptr),
    bandWidth_(dict.getOrDefault<label>("RDFbandWidth", 4)),
    curvSmoothCoeff_(dict.getOrDefault<scalar>("RDFcurvSmoothCoeff", 0.5)),
    nSmoothIter_(dict.getOrDefault<label>("RDFnSmoothIter", 2)),
    sigma_("sigma", dimensionSet(1, 0, -2, 0, 0, 0, 0), 0.0)
{
    // Read surface tension from constant/transportProperties. We reach it
    // via a fresh IOdictionary handle rather than mixture_.sigma() so we
    // do not depend on whether interfaceProperties is publicly reachable
    // on the mixture. Two layouts are supported:
    //   (a) top-level:            sigma   0.072;
    //   (b) sub-dict:             sigma   sigma [1 0 -2 0 0] 0.072;
    //       (the format immiscibleIncompressibleTwoPhaseMixture itself
    //        uses internally when parsing transportProperties)
    {
        IOdictionary transportDict
        (
            IOobject
            (
                "transportProperties",
                alpha1.mesh().time().constant(),
                alpha1.mesh(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );

if (transportDict.found("sigma"))
{
    // v2212 现代 API:第三个参数是显式量纲,从 dict 读 key="sigma".
    // 同时兼容纯标量 "sigma 0.072;" 和完整形式 "sigma sigma [1 0 -2 0 0 0 0] 0.072;"
    sigma_ = dimensionedScalar
    (
        "sigma",
        dimensionSet(1, 0, -2, 0, 0, 0, 0),
        transportDict
    );
}
else
{
    FatalErrorInFunction
        << "Cannot find 'sigma' entry in constant/transportProperties."
        << exit(FatalError);
}
    }

    Info<< "Constructing RDF curvature helper" << nl
        << "    sigma               = " << sigma_.value() << " N/m" << nl
        << "    RDFbandWidth        = " << bandWidth_       << nl
        << "    RDFcurvSmoothCoeff  = " << curvSmoothCoeff_ << nl
        << "    RDFnSmoothIter      = " << nSmoothIter_     << nl
        << endl;

    // The isoAdvection constructor reads "alpha.<phaseName>" sub-dict
    // from system/fvSolution (same as interIsoFoam). For a case where
    // only MPPICInterFoam-style settings exist, isoAdvection will fall
    // back to its defaults, but it is good practice to include an
    // isoAdvector sub-dict. See README.
    advectorPtr_.reset(new isoAdvection(alpha1_, phi_, U_));

    // Construct the RDF object (size = nCells, values filled per update)
    rdfPtr_.reset(new reconstructedDistanceFunction(mesh_));

    // Initial reconstruction so kappa_RDF is well-defined on the first
    // UEqn assembly, before alphaEqnSubCycle.H has run.
    update();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigmaKRDF::update()
{
    // 1. Geometric iso-surface reconstruction (PLIC + RDF by default).
    //    This populates surf().interfaceCell(), centre(), normal().
    advectorPtr_->surf().reconstruct();

    // 2. Mark cells within bandWidth_ neighbours of an interface cell.
    rdfPtr_->markCellsNearSurf
    (
        advectorPtr_->surf().interfaceCell(),
        bandWidth_
    );

    // 3. Build the signed distance field via v2212 API.
    //    zoneDistribute is a registered singleton per mesh.
    zoneDistribute& exchange = zoneDistribute::New(mesh_);
    exchange.setUpCommforZone(rdfPtr_->nextToInterface(), true);

    rdfPtr_->constructRDF
    (
        rdfPtr_->nextToInterface(),
        advectorPtr_->surf().centre(),
        advectorPtr_->surf().normal(),
        exchange
    );

    // Keep boundary values physically reasonable. RDF extrapolates into
    // boundary cells but we do not rely on those values: kappa_RDF is
    // only trusted in the narrow-band mask.
    rdfPtr_->correctBoundaryConditions();
}


Foam::tmp<Foam::volScalarField> Foam::sigmaKRDF::sigmaK() const
{
    const reconstructedDistanceFunction& rdf = rdfPtr_();

    // Small deltaN, same recipe used by interfaceProperties.
const dimensionedScalar deltaN
(
    "deltaN",
    dimless,
    1e-8
);

    // -------- RDF-based curvature ---------
    //
    // kappa_RDF = -div( grad(RDF) / |grad(RDF)| )
    //
    // Grad(RDF) should be ~unit vector in the narrow band (since RDF is
    // a signed distance). Normalising it explicitly guards against tiny
    // numerical drift.
    tmp<volVectorField> tgRDF = fvc::grad(rdf);
    const volVectorField& gRDF = tgRDF();

    volVectorField nRDF
    (
        IOobject
        (
            "nRDF",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        gRDF/(mag(gRDF) + deltaN)
    );
    nRDF.correctBoundaryConditions();

    volScalarField kRDF
    (
        IOobject
        (
            "kappaRDF",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        -fvc::div(nRDF)
    );

    // Optional Laplacian smoothing of kappa_RDF to suppress residual
    // high-frequency noise when the interface straddles cells of very
    // different aspect ratio.
    for (label i = 0; i < nSmoothIter_; ++i)
    {
        kRDF =
            (1.0 - curvSmoothCoeff_)*kRDF
          + curvSmoothCoeff_*fvc::average(fvc::interpolate(kRDF));
    }

    // -------- Narrow-band mask ---------
    //
    // The reconstructedDistanceFunction exposes nextToInterface() which
    // is true in all cells that received a valid distance. We build a
    // scalar mask from it and soften the transition with one Laplacian
    // smoothing pass so the blend does not introduce a visible step in
    // sigmaK.
    const boolList& inBand = rdf.nextToInterface();

    volScalarField mask
    (
        IOobject
        (
            "kappaMask",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );

    scalarField& m = mask.primitiveFieldRef();
    forAll(m, celli)
    {
        m[celli] = inBand[celli] ? 1.0 : 0.0;
    }
    mask.correctBoundaryConditions();

    // One pass of Laplacian smoothing to avoid a discontinuous blend
    mask = fvc::average(fvc::interpolate(mask));
    mask = min(max(mask, scalar(0)), scalar(1));

    // -------- Blend at the (sigma * kappa) level ---------
    //
    //  sigmaK_blended  =  mask * (sigma_ * kappa_RDF)
    //                  + (1 - mask) * mixture.sigmaK()
    //
    // Notes on why we blend at this level rather than at kappa:
    //   * mixture.sigmaK() is the ONLY public accessor that survives
    //     across OpenFOAM versions (some versions hide K() via private
    //     inheritance of interfaceProperties). Using it directly here
    //     avoids relying on mixture.K() / mixture.sigma() being
    //     publicly reachable.
    //   * sigma_ is read from transportProperties in the constructor
    //     so that we can multiply kappa_RDF into a field with the
    //     correct dimensions [1 -1 -2 0 0 0 0].
    //
    return tmp<volScalarField>::New
    (
        IOobject
        (
            "sigmaKRDF",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mask*(sigma_*kRDF) + (1.0 - mask)*mixture_.sigmaK()
    );
}


// ************************************************************************* //
