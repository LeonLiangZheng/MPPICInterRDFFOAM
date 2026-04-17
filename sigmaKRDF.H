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

Class
    Foam::sigmaKRDF

Description
    Helper that computes a surface-tension coefficient field

        sigmaK  =  sigma * kappa

    using curvature extracted from a reconstructed distance function (RDF),
    as an alternative to the algebraic CSF form used by
    interfaceProperties::sigmaK():

        kappa_CSF = -div( grad(alpha) / |grad(alpha)| )

    Strategy
    --------
    An auxiliary isoAdvection object is held. Its advect() method is
    NEVER called: alpha transport continues to be done by MULES in the
    outer solver so that the MPPIC-cloud coupling (alphac, alphaSuSp,
    Sp + divU) is unchanged.

    On each call to update(), we:

      1. Run surf().reconstruct() to refresh the iso-surface from the
         current alpha1 (built by MULES in the preceding step).
      2. Ask the reconstructionSchemes for its interfaceCell() list,
         centre() and normal() fields.
      3. Build a reconstructedDistanceFunction (RDF) in a narrow band
         around the interface using the v2212 API:

            rdf.markCellsNearSurf(interfaceCell_, bandWidth);
            rdf.constructRDF(nextToInterface, centre, normal, exchange);

      4. Compute kappa_RDF = -div( grad(RDF) / |grad(RDF)| ) and
         optionally smooth it by Laplacian averaging.

    On each call to sigmaK(), we build a blended curvature

        kappa = m * kappa_RDF + (1 - m) * kappa_CSF

    where m is a narrow-band mask based on reconstructedDistanceFunction
    membership (rdf.nextToInterface()). Outside the band we fall back
    onto the algebraic CSF curvature (which is already stored in
    mixture.K()), so the contact-angle BC carried by mixture.nHatf() is
    preserved at walls intersecting the interface.

    Output has the same dimensions as mixture.sigmaK() and can be used
    as a drop-in replacement in UEqn.H and pEqn.H.

Caveats
    - plicRDF is the default isoAdvection reconstruction scheme in
      v2212 ESI. It exposes centre() and normal() as volVectorFields on
      interfaceCells (zero elsewhere), which is exactly what
      reconstructedDistanceFunction::constructRDF() expects.

    - zoneDistribute is a singleton-style parallel stencil cache. We
      fetch it via zoneDistribute::New(mesh) and re-configure its zone
      before each RDF construction.

    - reconstructedDistanceFunction is itself a volScalarField (it
      derives from volScalarField). Values OUTSIDE the narrow band are
      arbitrary; do not use them.

SourceFiles
    sigmaKRDF.C

\*---------------------------------------------------------------------------*/

#ifndef sigmaKRDF_H
#define sigmaKRDF_H

#include "fvCFD.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "isoAdvection.H"
#include "reconstructedDistanceFunction.H"
#include "zoneDistribute.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class sigmaKRDF Declaration
\*---------------------------------------------------------------------------*/

class sigmaKRDF
{
    // Private data

        //- Mesh
        const fvMesh& mesh_;

        //- Mixture (for sigma and for the CSF fallback curvature)
        immiscibleIncompressibleTwoPhaseMixture& mixture_;

        //- alpha1 (passed through to isoAdvection)
        volScalarField& alpha1_;

        //- Flux and velocity (passed through to isoAdvection)
        const surfaceScalarField& phi_;
        const volVectorField& U_;

        //- Auxiliary isoAdvection. Used only for geometric reconstruction.
        autoPtr<isoAdvection> advectorPtr_;

        //- The signed distance field, owned here so it persists between
        //  reconstructions.
        autoPtr<reconstructedDistanceFunction> rdfPtr_;

        //- Narrow-band half-width (in cells) for RDF validity
        label bandWidth_;

        //- Laplacian smoothing of kappa_RDF
        scalar curvSmoothCoeff_;
        label  nSmoothIter_;

        //- Surface tension. Read from constant/transportProperties in the
        //  constructor. We store it here because mixture_.sigma() is not
        //  reliably publicly reachable on immiscibleIncompressibleTwoPhaseMixture
        //  across OpenFOAM versions.
        dimensionedScalar sigma_;


public:

    // Constructor
        sigmaKRDF
        (
            immiscibleIncompressibleTwoPhaseMixture& mixture,
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U,
            const dictionary& dict
        );


    // Member functions

        //- Refresh the iso-surface reconstruction and the RDF narrow band.
        //  Call this after MULES has produced an updated alpha1.
        void update();

        //- Return sigma * kappa on cell centres, blended between
        //  RDF-based (near interface) and CSF-based (elsewhere) curvature.
        tmp<volScalarField> sigmaK() const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
