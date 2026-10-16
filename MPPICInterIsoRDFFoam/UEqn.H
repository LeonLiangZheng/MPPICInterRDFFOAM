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
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            phicForces/rAUcf
          +
            (
                fvc::interpolate
                (
                    // 【方案 B 修改点】曲率源头：RDF 模式 → sigmaKRDF；CSF 模式 → mixture.sigmaK()
                    useRDFCurvature
                        ? rdfCurvPtr->sigmaK()
                        : mixture.sigmaK()
                )*fvc::snGrad(alpha1)
                - ghf*fvc::snGrad(rho)
                - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
    );

    fvOptions.correct(U);
}
