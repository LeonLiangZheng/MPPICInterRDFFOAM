/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | www.openfoam.com
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Derived from OpenFOAM v2212 (ESI) solvers:
      - MPPICInterFoam  (Copyright (C) 2016-2020 OpenCFD Ltd.)
      - interIsoFoam    (Copyright (C) 2011-2020 OpenFOAM Foundation / OpenCFD /
                         DHI / Johan Roenby / DLR)
-------------------------------------------------------------------------------
License
    GPL-3.0-or-later (inherits from upstream OpenFOAM components).

Application
    MPPICInterRDFFoam

Description
    Two-phase incompressible VOF solver with an MPPIC Lagrangian cloud
    (as in MPPICInterFoam), extended with Reconstructed-Distance-Function
    (RDF) interface-curvature evaluation borrowed from interIsoFoam's
    geometricVoF framework.

    Rationale
    ---------
    Algebraic CSF curvature kappa = -div(grad(alpha)/|grad(alpha)|) is a
    second-order numerical derivative of a near-discontinuous field and
    produces large spurious currents near contact lines in droplet-on-wire
    type problems (Weber ~ O(10^2)-O(10^3)).

    This solver keeps the MULES-based alpha transport of MPPICInterFoam
    (so that the alphac <-> kinematicCloud coupling is untouched) and
    adds an auxiliary isoAdvection object whose geometric interface
    reconstruction yields a smooth signed-distance field (RDF) in a narrow
    band around the interface. Curvature is then computed from the RDF,
    which is smoother than alpha by construction.

    At run-time, select the curvature model in
        system/fvSolution -> PIMPLE -> curvatureModel
    Supported values:
        CSF  (default)   : original algebraic CSF, identical to
                           MPPICInterFoam
        RDF              : RDF-based curvature in the narrow band,
                           blended with CSF outside

    With RDF selected, the surface-tension force used in UEqn and pEqn is
        F_sigma = sigma * kappa_RDF * grad(alpha)
    i.e. the CSF *form* is preserved (alpha is still the indicator used
    for force application), but kappa comes from the smooth RDF.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"

#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"

#include "basicKinematicCloud.H"

// --- geometricVoF / RDF
#include "isoAdvection.H"
#include "reconstructedDistanceFunction.H"

// --- Local helper for RDF-based curvature
#include "sigmaKRDF.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Two-phase VOF solver (MULES) with MPPIC Lagrangian cloud,"
        " extended with RDF-based interface curvature from geometricVoF."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createAlphaFluxes.H"
    #include "createFvOptions.H"
    #include "correctPhi.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        Info<< "Evolving " << kinematicCloud.name() << endl;

        kinematicCloud.evolve();

        // Update continuous phase volume fraction field
        alphac = max(1.0 - kinematicCloud.theta(), alphacMin);
        alphac.correctBoundaryConditions();

        Info<< "Continuous phase-1 volume fraction = "
            << alphac.weightedAverage(mesh.Vsc()).value()
            << "  Min(alphac) = " << min(alphac).value()
            << "  Max(alphac) = " << max(alphac).value()
            << endl;

        alphacf = fvc::interpolate(alphac);
        alphaRhoPhic = alphacf*rhoPhi;
        alphaPhic = alphacf*phi;
        alphacRho = alphac*rho;

        fvVectorMatrix cloudSU(kinematicCloud.SU(U));
        volVectorField cloudVolSUSu
        (
            IOobject
            (
                "cloudVolSUSu",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedVector(cloudSU.dimensions()/dimVolume, Zero),
            zeroGradientFvPatchVectorField::typeName
        );

        cloudVolSUSu.primitiveFieldRef() = -cloudSU.source()/mesh.V();
        cloudVolSUSu.correctBoundaryConditions();

        cloudSU.source() = vector::zero;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

            // After MULES has updated alpha1, refresh the RDF narrow band.
            // We only call this when RDF curvature is actually requested,
            // as the geometric reconstruction is not free.
            if (useRDFCurvature)
            {
                rdfCurvPtr->update();
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
