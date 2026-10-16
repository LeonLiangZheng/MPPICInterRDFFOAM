# MPPICInterRDFFoam (legacy)

> **Superseded.** This is the first (April 2026) version of the solver, kept for
> reference. It keeps MULES transport and runs an auxiliary isoAdvection object only
> to reconstruct the interface. The current solvers are in
> [`MPPICInterIsoRDFFoam`](../../MPPICInterIsoRDFFoam) and
> [`MPPICInterIsoFoam`](../../MPPICInterIsoFoam); see the [top-level README](../../README.md).

MPPICInterFoam with RDF-based interface curvature, derived from two OpenFOAM
v2212 (ESI) solvers:

- `MPPICInterFoam`: MULES VOF + MPPIC Lagrangian cloud
- `interIsoFoam`: isoAdvection geometric VOF + reconstructedDistanceFunction

The alpha transport in this solver remains MULES (identical to MPPICInterFoam)
so that the continuous-phase fraction coupling with the kinematicCloud is
untouched. What this solver adds is an auxiliary isoAdvection object used
purely for geometric interface reconstruction, from which a smooth
reconstructed distance function (RDF) is obtained. The interface curvature
used in the CSF surface-tension source term is then computed as

    kappa_RDF  =  -div( grad(RDF) / |grad(RDF)| )

in cells that are inside the RDF narrow band, and falls back to the
algebraic CSF curvature

    kappa_CSF  =  -div( grad(alpha) / |grad(alpha)| )      ( = mixture.K() )

elsewhere. This preserves the Kistler dynamic contact-angle treatment at
wall-adjacent cells, which is carried by `mixture.nHatf()` in CSF.

## 1. Directory layout

```
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
`-- compressibleTwoPhaseMixtureTurbulenceModels/
                                   (the same library MPPICInterFoam ships)
```

## 2. Build

```bash
source /path/to/OpenFOAM-v2212/etc/bashrc
cd legacy/MPPICInterRDFFoam
./Allwmake
```

or, equivalently:

```bash
wmake libso compressibleTwoPhaseMixtureTurbulenceModels
wmake
```

The first line builds the mixture-turbulence library (the same library that
MPPICInterFoam ships with; if MPPICInterFoam is already installed, its
`libcompressibleTwoPhaseMixtureTurbulenceModels` is reused and this step is a
no-op). The second line builds the solver itself into `$FOAM_USER_APPBIN`.
To clean, run `./Allwclean`.

## 3. Required case-file changes

### 3.1 `system/fvSolution`, PIMPLE

Add the curvature-model keyword and (if using the RDF path) the `*Final`
entries required by v2212 when `nCorrectors > 1`.

```
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
```

When `curvatureModel` is absent or set to `CSF`, the solver is numerically
equivalent to MPPICInterFoam; the RDF machinery is not even instantiated.

### 3.2 `system/fvSolution`, alpha solver block

The auxiliary isoAdvection object reads its controls from the alpha solver
block, in the same way interIsoFoam does. Even though isoAdvector is not used
for transport, include these entries for robust reconstruction:

```
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
```

### 3.3 `system/fvSchemes`

No changes required beyond what MPPICInterFoam already needs. The extra
`fvc::grad(rdf)` and `fvc::div(nRDF)` operators in `sigmaKRDF.C` pick up the
default schemes; to tune them, add:

```
gradSchemes
{
    grad(reconstructedDistanceFunction) Gauss linear;
    grad(nRDF)                           Gauss linear;
}

divSchemes
{
    div(nRDF) Gauss linear;
}
```

## 4. Recommended first tests

1. **Baseline**: build, set `curvatureModel CSF`, run an existing
   MPPICInterFoam case unchanged. Output should match MPPICInterFoam up to
   round-off; if not, there is a build or include issue before RDF is even on
   the table.
2. **Static droplet**: a sphere of radius R in zero gravity, no particles, zero
   initial velocity, `curvatureModel RDF`. Inspect `kappaRDF`, `kappaBlend` and
   `sigmaKRDF`: the droplet should show kappa ≈ 2/R with a spread of a few
   percent (CSF typically scatters by 20–40 % on the same mesh), and spurious
   currents should drop by one to two orders of magnitude.
3. **Short production test**: the real case with `curvatureModel RDF` on a
   coarse mesh for about 10 time steps; compare the peak `p_rgh` against the
   CSF run.
4. Only then deploy to the cluster.

## 5. Known caveats

- RDF validity is limited to a narrow band. Cells outside it have undefined
  RDF values; the blending mask makes those cells use CSF curvature. The band
  is controlled by `RDFbandWidth`.
- After strong topology changes (pinch-off, satellite droplets) the
  reconstruction can produce multi-valued RDF near a just-formed daughter
  droplet, and kappa_RDF may briefly be noisy. The Laplacian smoothing pass
  (`RDFcurvSmoothCoeff`, `RDFnSmoothIter`) is the remedy; if noise persists,
  increase `RDFnSmoothIter` to 3 or 4.
- Where the interface meets a wall patch, `grad(RDF)` is ill-defined. The
  Kistler contact-angle BC acts on alpha through `mixture.nHatf()` in the CSF
  branch; because the band mask excludes many wall-adjacent cells, the blended
  curvature uses CSF there and the contact-angle treatment is preserved. If
  the contact line drifts, reduce `RDFbandWidth` (e.g. to 2).
- Parallel: `zoneDistribute` and `reconstructedDistanceFunction` are
  parallel-aware in v2212, but the solver has not been validated on more than
  8 ranks; start a large run with a short (10-step) test to rule out deadlocks
  in `zoneDistribute::setUpCommforZone`.
- No dynamic mesh refinement: that would need `sigmaKRDF::update()` to rebuild
  the RDF stencil after topology changes.

## 6. Dimensional check

```
sigma              [1 0 -2 0 0 0 0]   (kg/s^2, i.e. N/m)
kappa              [0 -1 0 0 0 0 0]   (1/length)
sigmaK = sigma*K   [1 -1 -2 0 0 0 0]
fvc::interpolate(sigmaK) * fvc::snGrad(alpha) * mesh.magSf()
                   [1 1 -2 0 0 0 0]   (N, a force)
```

This matches the form used by `mixture.surfaceTensionForce()` in stock
interFoam / interIsoFoam, so no scaling fix-ups are needed.

## 7. If the build fails

Most likely causes, in order of frequency:

- **Missing geometricVoF library.** Check
  `ls $FOAM_SRC/transportModels/geometricVoF/lnInclude/isoAdvection.H`. It
  ships with v2212 ESI by default but may be stripped from a minimal install.
- **Linker error on `-lgeometricVoF`.** Confirm
  `ls $FOAM_LIBBIN/libgeometricVoF.so`; if present, try moving
  `-lgeometricVoF` earlier in `EXE_LIBS`.
- **"unresolved reference to zoneDistribute::New".** The class is in
  libgeometricVoF in v2212; on another version, find it with
  `grep -r "class zoneDistribute" $FOAM_SRC/` and update `Make/options`.
- **Warning about an inheritance diamond for reconstructedDistanceFunction.**
  Expected (the class derives publicly from volScalarField); no action needed
  unless `-Werror` is set.
