/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

Application
    MPPICInterIsoFoam

Description
    方案 A 衍生求解器：在 MPPICInterFoam 基础上，把 MULES 界面输运替换为
    isoAdvector 几何输运。动量方程、压力方程、MPPIC 粒子云、曲率（CSF）
    全部保留原 MPPICInterFoam 行为。

    核心修改点集中在 alphaEqn.H：对含 αc 的原 α 方程
        ∂(αc·α1)/∂t + ∇·(αc·α1·U) = Su + Sp·α1
    作链式法则变形得到
        ∂α1/∂t + ∇·(α1·U) = Su/αc − (α1/αc)·[∂αc/∂t + ∇·(αc·U)]
    使得 isoAdvector::advect(Sp, Su) 可直接吸收变形后的源项。

    保留原 MPPICInterFoam 的：
      - αc 加权 UEqn.H (fvm::ddt(alphacRho, U) + fvm::div(rhoPhi, U))
      - αc 加权 pEqn.H (Dp = alphacf*rAUcf, phi = phiHbyA − p_rghEqn.flux()/alphacf)
      - MPPIC cloud evolve 在 PIMPLE 循环外
      - kOmegaSST / PhaseCompressibleTurbulenceModel 框架
      - mixture.sigmaK() CSF 曲率（本方案不动）

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
// 【方案 A 新增】isoAdvector 头文件
#include "isoAdvection.H"
// 保留原 MPPICInterFoam 的 CMULES / CrankNicolson 头文件以保证兼容头内的
// ddtScheme 选择工作正常（alphaSuSp.H 等会间接用到）
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "MPPICInterFoam 衍生求解器：以 isoAdvector 几何输运替换 MULES；"
        "保留 MPPIC cloud、αc 加权动量/压力方程、CSF 曲率。"
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
    // createAlphaFluxes.H 仍然创建 alphaPhi10（alphaCourantNo.H 要读）
    // 但本求解器的 alphaPhi10 将在 alphaEqn.H 里被赋值为 αcf*advector.alphaPhi()
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
            // 【注意】alphaCourantNo.H 用的是 alphaPhi10，在方案 A 下
            // alphaPhi10 由 αcf*advector.alphaPhi() 提供（见 alphaEqn.H）。
            // 建议 case 中设置 maxAlphaCo ≈ 0.3（isoAdvector 比 MULES 严：
            // 稳定极限约 0.5）
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        Info<< "Evolving " << kinematicCloud.name() << endl;

        // MPPIC 粒子云在 PIMPLE 之外 evolve —— 进入 PIMPLE 时 αc 已固定
        kinematicCloud.evolve();

        // 更新连续相体积分数场
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
