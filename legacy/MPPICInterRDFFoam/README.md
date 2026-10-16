> **Superseded.** This is the first (April 2026) version of the solver, kept for
> reference. It keeps MULES transport and runs an auxiliary isoAdvection object only
> to reconstruct the interface. The current solvers are in
> [`../../MPPICInterIsoRDFFoam`](../../MPPICInterIsoRDFFoam) and
> [`../../MPPICInterIsoFoam`](../../MPPICInterIsoFoam); see the top-level README.

================================================================================
  MPPICInterRDFFoam  —  MPPICInterFoam with RDF-based interface curvature
================================================================================

Derived from two OpenFOAM v2212 (ESI) solvers:
  - MPPICInterFoam  : MULES VOF + MPPIC Lagrangian cloud
  - interIsoFoam    : isoAdvection geometric VOF + reconstructedDistanceFunction

The alpha transport in this solver remains MULES (identical to MPPICInterFoam)
so that the continuous-phase fraction coupling with the kinematicCloud is
untouched. What this solver ADDS is an auxiliary isoAdvection object used
purely for geometric interface reconstruction, from which a smooth
reconstructed distance function (RDF) is obtained. The interface curvature
kappa used in the CSF surface-tension source term is then computed as

    kappa_RDF  =  -div( grad(RDF) / |grad(RDF)| )

in cells that are inside the RDF narrow band, and falls back to the
algebraic CSF curvature

    kappa_CSF  =  -div( grad(alpha) / |grad(alpha)| )      ( = mixture.K() )

elsewhere. This preserves the Kistler dynamic contact-angle treatment at
wall-adjacent cells, which is carried by mixture.nHatf() in CSF.


--------------------------------------------------------------------------------
1.  Directory layout
--------------------------------------------------------------------------------
MPPICInterRDFFoam/
|-- MPPICInterRDFFoam.C            main solver
|-- sigmaKRDF.H                    RDF curvature helper (declaration)
|-- sigmaKRDF.C                    RDF curvature helper (implementation)
|-- createFields.H                 fields + reads curvatureModel from fvSolution
|-- UEqn.H, pEqn.H                 momentum / pressure equations (curv switch)
|-- alphaEqn.H, alphaEqnSubCycle.H unchanged MULES path (verbatim from MPPIC)
|-- continuityErrs.H, correctPhi.H unchanged
|-- Make/files, Make/options       build rules
|-- Allwmake, Allwclean            convenience scripts
|-- compressibleTwoPhaseMixtureTurbulenceModels/
|   `-- ...                        (the same library MPPICInterFoam ships)


--------------------------------------------------------------------------------
2.  Build
--------------------------------------------------------------------------------
Source OpenFOAM-v2212 first:

    source /path/to/OpenFOAM-v2212/etc/bashrc

then from inside the MPPICInterRDFFoam/ directory:

    ./Allwmake

or, equivalently:

    wmake libso compressibleTwoPhaseMixtureTurbulenceModels
    wmake

The first line builds the mixture-turbulence library (this is the same
library that MPPICInterFoam ships with; if MPPICInterFoam is already
installed on the system, its libcompressibleTwoPhaseMixtureTurbulenceModels
is reused and this step is a no-op). The second line builds the solver
itself into $FOAM_USER_APPBIN.

To clean:

    ./Allwclean


--------------------------------------------------------------------------------
3.  Required case-file changes
--------------------------------------------------------------------------------

3.1  system/fvSolution
------------------------

Add two things inside the PIMPLE sub-dictionary: the curvature-model
keyword and (if using the RDF path) the *Final entries required by v2212
when nCorrectors > 1.

    PIMPLE
    {
        momentumPredictor            yes;
        nOuterCorrectors             3;
        nCorrectors                  3;
        nNonOrthogonalCorrectors     0;

        // --- New keywords ---
        curvatureModel               RDF;   // or CSF (default = original behaviour)
        RDFbandWidth                 4;     // narrow band, in cell layers (default 4)
        RDFcurvSmoothCoeff           0.5;   // 0..1 Laplacian smoothing (default 0.5)
        RDFnSmoothIter               2;     // smoothing passes on kappa_RDF (default 2)
    }

When curvatureModel is absent or set to CSF, the solver is byte-for-byte
equivalent in numerics to MPPICInterFoam — the RDF machinery is not
even instantiated.


3.2  system/fvSolution (isoAdvector sub-dict)
----------------------------------------------

The auxiliary isoAdvection object reads an "isoAdvector" sub-dict inside
the alpha1 solver block, in the same way interIsoFoam does. Even if you
are not using isoAdvector for transport, include this block for robust
reconstruction:

    "alpha.water.*"
    {
        // --- existing MULES entries (keep them) ---
        nAlphaCorr       2;
        nAlphaSubCycles  2;
        cAlpha           2;
        MULESCorr        no;
        nLimiterIter     10;
        solver           smoothSolver;
        smoother         symGaussSeidel;
        tolerance        1e-8;
        relTol           0;
        maxIter          100;

        // --- isoAdvector reconstruction controls ---
        isoFaceTol       1e-10;
        surfCellTol      1e-8;
        nAlphaBounds     3;
        snapTol          1e-12;
        clip             true;

        // --- PLIC reconstruction scheme (default is plicRDF) ---
        reconstructionScheme plicRDF;  // or isoAlpha, gradAlpha
    }


3.3  system/fvSchemes
----------------------

No changes required beyond what MPPICInterFoam already needs. The extra
fvc::grad(rdf) and fvc::div(nRDF) operators in sigmaKRDF.C pick up the
default ddt / grad / div schemes; if you want to tune them, add:

    gradSchemes
    {
        grad(reconstructedDistanceFunction) Gauss linear;
        grad(nRDF)                           Gauss linear;
    }

    divSchemes
    {
        div(nRDF) Gauss linear;
    }


--------------------------------------------------------------------------------
4.  Recommended first tests (on the VM, before going to HPC)
--------------------------------------------------------------------------------

Step 0   Baseline sanity — build with curvatureModel = CSF, run your
         existing MPPICInterFoam case unchanged. Output should be
         byte-identical (up to round-off) to MPPICInterFoam output. If
         this step disagrees, there is a build/include issue before RDF
         is even on the table.

Step 1   Static droplet test — put a sphere of radius R in zero gravity,
         no MPPIC particles, zero initial velocity. Run with:
             curvatureModel  RDF;
         and inspect kappaRDF, kappaBlend, sigmaKRDF in the output. The
         interior of the droplet should show kappa ~ 2/R with spread of
         a few percent; in MPPICInterFoam (CSF) the same field typically
         shows 20-40 % scatter on the same mesh. Spurious currents should
         drop by 1-2 orders of magnitude.

Step 2   Your actual droplet-on-wire case, curvatureModel = RDF, on the
         VM with the existing coarse mesh, ~10 time steps. Compare the
         p_rgh peak magnitude at t ~ 0.5 ms against the CSF run:
             CSF:  ~23800 Pa  (previous observation)
             RDF:  target <~ 2000 Pa, well below MPPIC-induced peaks

Step 3   Only then deploy to the HPC cluster with 64 ranks.


--------------------------------------------------------------------------------
5.  Known caveats
--------------------------------------------------------------------------------

- RDF validity is limited to a narrow band. Cells outside that band have
  undefined RDF values; the blending mask ensures those cells use CSF
  curvature instead. The band is controlled by RDFbandWidth.

- After strong topology changes (droplet pinch-off, satellite formation)
  the iso-surface reconstruction can produce multi-valued RDF in cells
  close to a just-formed daughter droplet. In such cells kappa_RDF may
  briefly be noisy. The Laplacian smoothing pass
  (RDFcurvSmoothCoeff, RDFnSmoothIter) is the remedy; if noise persists,
  increase RDFnSmoothIter to 3 or 4.

- In cells where the interface coincides with a wall patch intersecting
  the geometry (the wire surface), grad(RDF) is ill-defined. The Kistler
  contact-angle BC is applied to alpha1 through mixture.nHatf() in the
  CSF branch. Because the RDF narrow-band mask naturally excludes many
  of these wall-adjacent cells, the blended kappa uses CSF there and
  the contact-angle treatment is preserved. If the contact line shows
  drift in validation, reduce RDFbandWidth (e.g. to 2) so that the mask
  becomes tighter.

- Parallel: zoneDistribute and reconstructedDistanceFunction are
  parallel-aware in v2212. The solver has not been validated on >8
  ranks by the author; the first production run on 64 ranks should be
  a short (10-step) test to confirm absence of deadlocks in
  zoneDistribute::setUpCommforZone for the specific decomposition.

- This solver does NOT support dynamic mesh refinement. Adding it
  would require additional logic in sigmaKRDF::update() to re-init the
  RDF stencil after topology changes.


--------------------------------------------------------------------------------
6.  Dimensional sanity check
--------------------------------------------------------------------------------

mixture.sigma()    has dimensions [1 0 -2 0 0 0 0]  (kg/s^2, i.e. N/m)
kappa              has dimensions [0 -1 0 0 0 0 0]  (1/length)
sigmaK = sigma*K   has dimensions [1 -1 -2 0 0 0 0]
fvc::interpolate(sigmaK) * fvc::snGrad(alpha) * mesh.magSf()
                   has dimensions [1 1 -2 0 0 0 0]  (N -- a force)

This matches exactly the form used by mixture.surfaceTensionForce() in
stock interFoam/interIsoFoam, so no scaling fix-ups are needed.


--------------------------------------------------------------------------------
7.  If the build fails
--------------------------------------------------------------------------------

Most likely causes, in order of frequency:

(a) Missing geometricVoF library. Check:
        ls $FOAM_SRC/transportModels/geometricVoF/lnInclude/isoAdvection.H
    If absent, your OpenFOAM install does not have the geometric VOF
    module. It ships with v2212 ESI by default, but may have been
    stripped from a minimal install.

(b) Linker error on -lgeometricVoF. Confirm:
        ls $FOAM_LIBBIN/libgeometricVoF.so
    If present but linker still fails, the ordering in Make/options
    may matter; try moving -lgeometricVoF earlier in EXE_LIBS.

(c) "unresolved reference to zoneDistribute::New". This class is in
    libgeometricVoF in v2212; if you are on an older v-series and
    zoneDistribute has been moved, check:
        grep -r "class zoneDistribute" $FOAM_SRC/
    and update the include path in Make/options accordingly.

(d) Compiler warning about inheritance diamond for
    reconstructedDistanceFunction. This is expected — the class derives
    publicly from volScalarField. No action required unless -Werror is
    set.


================================================================================
  End of README
================================================================================
