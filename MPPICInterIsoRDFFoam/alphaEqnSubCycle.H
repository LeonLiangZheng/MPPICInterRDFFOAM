// ============================================================================
// alphaEqnSubCycle.H （方案 A）
// ----------------------------------------------------------------------------
// 合并两个原版的核心逻辑：
//   1) interIsoFoam 的 prevIter 防震荡技巧（nOuterCorrectors > 1 时必要）
//   2) MPPICInterFoam 的 subCycle + rhoPhiSum 权重平均
//
// prevIter 技巧说明：
//   当 nOuterCorrectors > 1 时，外迭代会重复解 α 方程，若每次使用新 U/phi
//   直接覆盖 α1 会导致界面震荡。策略：
//     - 首次外迭代：存储 α1、U、phi 的初值
//     - 后续外迭代：α1 回滚到初值；U、phi 换为与当前步一半平均
//   之后再复原 U、phi 供 UEqn/pEqn 使用。
//
//   当 nOuterCorrectors == 1（Leon 正在调的性能方向），
//   `if (pimple.firstIter())` 永远 true，else 分支从不执行，开销为零。
// ============================================================================

if (pimple.nCorrPIMPLE() > 1)
{
    if (pimple.firstIter())
    {
        alpha1.storePrevIter();
        phi.storePrevIter();
        U.storePrevIter();
    }
    else
    {
        alpha1 = alpha1.prevIter();
        U = 0.5*U.prevIter() + 0.5*U;
        phi = 0.5*phi.prevIter() + 0.5*phi;
    }
}

if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(rhoPhi.dimensions(), Zero)
    );

    // 注：LTS 下 isoAdvector 的子步 deltaT 由 subCycle 自行拆分；不像 MULES
    // 那样需要显式传 trSubDeltaT。保留 LTS 兼容性无需额外处理。

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
    #include "alphaEqn.H"
}

// prevIter 恢复 U、phi 到本迭代真值（供 UEqn/pEqn 使用）
if (pimple.nCorrPIMPLE() > 1 && !pimple.firstIter())
{
    U = 2.0*U - U.prevIter();
    phi = 2.0*phi - phi.prevIter();
}

rho == alpha1*rho1 + alpha2*rho2;
mu = mixture.mu();
