# MPPICInterRDFFoam

OpenFOAM v2212 (ESI) solvers for two-phase VOF flow that carries an MPPIC
particle cloud. They start from the stock `MPPICInterFoam` and change two things:

- **Interface transport**: MULES is replaced by isoAdvector geometric transport,
  so the interface stays sharp instead of being smeared over several cells.
- **Surface-tension curvature** (optional): near the interface, curvature is
  computed from a reconstructed distance function (RDF) instead of the CSF
  formula, which reduces spurious velocities around the interface and contact
  lines.

## What is in this repository

| Folder | Status | What it is |
|---|---|---|
| [`MPPICInterIsoRDFFoam/`](MPPICInterIsoRDFFoam) | current | isoAdvector transport + curvature switchable between RDF and CSF at run time |
| [`MPPICInterIsoFoam/`](MPPICInterIsoFoam) | current | isoAdvector transport, CSF curvature only; the baseline for A/B comparison |
| [`legacy/MPPICInterRDFFoam/`](legacy/MPPICInterRDFFoam) | superseded | the first version (April 2026): MULES transport, RDF only for curvature |

Companion libraries, often used with these solvers:
[modifiedKistlerDynamicContactAngle](https://github.com/LeonLiangZheng/ModifiedKistlerDynamicContactAngle)
(dynamic contact-angle wall condition) and
[KOmegaSST4MPPICInterFOAM](https://github.com/LeonLiangZheng/KOmegaSST4MPPICInterFOAM)
(makes `kOmegaSST` selectable in the MPPICInterFoam family).

## How it works

**Transport.** `isoAdvection::advect(Sp, Su)` assumes every cell is fully
occupied by fluid, but with an MPPIC cloud the continuous phase only fills a
fraction `αc` of each cell, and the volume-fraction equation is

    ∂(αc·α)/∂t + ∇·(αc·α·U) = Su + Sp·α

Expanding with the chain rule (∇·U = 0) and dividing by `αc` gives an equation
isoAdvector can solve directly:

    ∂α/∂t + ∇·(α·U) = Su/αc + (Sp/αc − D/αc)·α,     D = ∂αc/∂t + ∇·(αc·U)

`MPPICInterIsoFoam` solves this equation as written; `MPPICInterIsoRDFFoam`
leaves out the `−D/αc` term, for the reason given [below](#the-d-term).
`αc` is clipped from below by `alphacMinClip` (default 0.05) wherever it is a
denominator. The mass flux is rebuilt as
`rhoPhi = αcf·alphaPhi·(ρ1 − ρ2) + αcf·phi·ρ2`, so it stays consistent with the
`αc`-weighted momentum and pressure equations of `MPPICInterFoam`.

**Curvature** (`MPPICInterIsoRDFFoam` only). After transport, the PLIC surface
that the advector has just reconstructed is reused (no second reconstruction)
to build an RDF in a narrow band around the interface. There

    κ_RDF = −∇·( ∇RDF / |∇RDF| )

and the surface-tension term is blended cell by cell,

    σκ = m·σ·κ_RDF + (1 − m)·(σκ)_CSF,     m = clamp((bandWidth − cellDistLevel) / bandWidth, 0, 1)

so interface cells use RDF and cells outside the band fall back to CSF. The
fallback keeps contact-angle boundary conditions working at walls. With
`curvatureModel CSF` (the default) the RDF object is never created.

## Build

```bash
source /path/to/OpenFOAM-v2212/etc/bashrc
cd MPPICInterIsoRDFFoam      # or MPPICInterIsoFoam
./Allwmake
```

`Allwmake` first builds the bundled `compressibleTwoPhaseMixtureTurbulenceModels`
library (the same one `MPPICInterFoam` ships), then the solver into
`$FOAM_USER_APPBIN`.

Do not add `-lincompressibleInterPhaseTransportModels` or
`-lVoFphaseTurbulentTransportModels` to `Make/options`: they register
turbulence models that clash with the bundled library and corrupt the heap at
exit.

## Case setup

`system/fvSolution`

```
solvers
{
    "alpha.water.*"
    {
        nAlphaSubCycles       2;
        reconstructionScheme  plicRDF;     // or isoAlpha
        vof2IsoTol            1e-8;
        surfCellTol           1e-6;
        nAlphaBounds          3;
        snapTol               1e-12;
        clip                  true;
        // alphacMinClip      0.05;        // optional
    }
}

PIMPLE
{
    // MPPICInterIsoRDFFoam only
    curvatureModel   RDF;       // RDF or CSF (default CSF)
    RDFsigma         0.0728;    // N/m, must equal sigma in transportProperties
    RDFbandWidth     4;         // band width in cell layers
    RDFnSmoothIter   0;         // smoothing of the RDF normal; keep 0 near contact lines
}
```

Remove the MULES keywords (`nAlphaCorr`, `cAlpha`, `MULESCorr`, `nLimiterIter`,
linear-solver entries) from the alpha block. isoAdvector silently ignores them,
which makes an old case look as if they still had an effect.

`system/fvSchemes`

```
divSchemes
{
    div(nRDF)    Gauss linear;     // MPPICInterIsoRDFFoam with RDF
    // div(phi,alpha) and div(phirb,alpha) are MULES entries and can go
}
gradSchemes
{
    grad(RDF)    Gauss linear;     // optional: unlimited gradient keeps the normal accurate
}
```

`system/controlDict`: keep `maxAlphaCo` at 0.3 or below. isoAdvector becomes
unstable above about 0.5, so the MULES habit of 1.0 does not carry over.

Keep surface tension consistent: `sigma` in `transportProperties`, `RDFsigma`,
and the `sigma` of a contact-angle condition such as
`modifiedKistlerAlphaContactAngle` must be the same value. Nothing checks this
at run time. `RDFsigma` exists because the v2212 mixture class does not expose
`sigma()`.

## Status

Verified:

- Both solvers compile with OpenFOAM v2212.
- `MPPICInterIsoFoam` runs a test case on an HPC cluster.
- `MPPICInterIsoRDFFoam` enters and advances the time loop on a multi-bubble
  reactor case, and has been run in RDF mode on a droplet-impact case, which is
  where the issue described under [The D term](#the-d-term) showed up.

Not yet verified:

- how much water mass drifts once the `D` term is left out;
- how much RDF reduces spurious currents compared with CSF (static-droplet test);
- RDF curvature when many bubbles split and merge: the PLIC surface jumps, and
  so does the RDF. If curvature looks unphysical, rerun with `curvatureModel CSF`
  to see whether RDF is the cause.

Other limits: RDF acts only inside the band; the small regularisation `deltaN`
is computed once from the mean cell size, so it is not updated under dynamic
mesh refinement.

## The D term

In RDF mode, the RDF field inside droplets showed faint speckles at the cells
that hold particles. The suspected cause is the term `D = ∂αc/∂t + ∇·(αc·U)`:
in exact arithmetic it is zero, but `αc` comes from projecting particles onto
the mesh rather than from a conservation equation, so `D` keeps a discrete
residual wherever particles cross cells. Fed into the α equation, it nudges α
away from 1 inside the liquid, PLIC then marks those cells as interface cells,
and the RDF is built around them.

Since 2026-05-19 `MPPICInterIsoRDFFoam` therefore leaves the `−D/αc` term out
of its `alphaEqn.H`. The price is that `αc·α` is no longer conserved exactly;
the drift should be small for dilute clouds but has not been measured, and
dense clouds (`αc` well below 1) may need a different fix. The two solvers'
`alphaEqn.H` files differ only in this term: to put it back, copy
`MPPICInterIsoFoam/alphaEqn.H` over `MPPICInterIsoRDFFoam/alphaEqn.H` and
rebuild.

## History

| Date | Change |
|---|---|
| 2026-04-16 | `MPPICInterRDFFoam`: MULES transport kept, an auxiliary isoAdvection object used only to reconstruct the interface for RDF curvature. Now in `legacy/`. |
| 2026-04-22 | Rewrite: `MPPICInterIsoFoam` (isoAdvector transport) and `MPPICInterIsoRDFFoam` (plus switchable RDF curvature sharing the same advector). |
| 2026-05-19 | `MPPICInterIsoRDFFoam` leaves the `−D/αc` term out of the α equation. |

## License

Derived from OpenFOAM's `MPPICInterFoam` and `interIsoFoam`, and therefore
GPL-3.0-or-later. Most source comments are in Chinese.
