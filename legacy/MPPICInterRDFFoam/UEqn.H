MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::ddt(alphacRho, U)
  + MRF.DDt(alphacRho, U)
  - fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), U)
  + fvm::div(rhoPhi, U)
  + turbulence->divDevRhoReff(U)
  ==
    fvOptions(rho, U)
  + cloudSU
);

UEqn.relax();
fvOptions.constrain(UEqn);

volScalarField rAUc(1.0/UEqn.A());
surfaceScalarField rAUcf(fvc::interpolate(rAUc));


surfaceScalarField phicForces
(
    (fvc::interpolate(rAUc*cloudVolSUSu) & mesh.Sf())
);


if (pimple.momentumPredictor())
{
    // Select curvature source: RDF-based (new) or algebraic CSF (original).
    // In both cases sigma*K is a cell-centred volScalarField; the surface
    // tension force in the CSF formulation is sigma*K * grad(alpha), so
    // the face contribution becomes fvc::interpolate(sigmaK)*snGrad(alpha).
    volScalarField sigmaKCurv
    (
        useRDFCurvature
      ? rdfCurvPtr->sigmaK()
      : mixture.sigmaK()
    );

    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            phicForces/rAUcf
          +
            (
                fvc::interpolate(sigmaKCurv)*fvc::snGrad(alpha1)
                - ghf*fvc::snGrad(rho)
                - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
    );

    fvOptions.correct(U);
}
