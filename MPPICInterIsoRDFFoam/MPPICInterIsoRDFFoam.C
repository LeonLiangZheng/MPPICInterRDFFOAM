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
    MPPICInterIsoRDFFoam

Description
    方案 B 衍生求解器：在 MPPICInterIsoFoam（方案 A）基础上增加 RDF 曲率
    可选路径。几何输运（isoAdvector + 方案 A αc 变形）、αc 加权 UEqn/pEqn、
    MPPIC 云 —— 全部继承方案 A 不变；唯一新增是 runtime 可切换的 RDF 曲率。

    curvatureModel（在 fvSolution 的 PIMPLE 子字典）：
      - RDF：窄带内 κ = -div(grad(RDF)/|grad(RDF)|)，窄带外 CSF fallback
      - CSF（默认）：全域 mixture.sigmaK()；此时 rdfCurvPtr 不构造，零开销

    架构（share 模式）：
      - sigmaKRDF 持主 advector 的 const 引用，复用 advect() 已做的 PLIC 重构
      - 不为 RDF 再跑一次 PLIC（避免 64 rank HPC 上白给开销）
      - update() 时序：alphaEqnSubCycle.H → rdfCurv.update() → UEqn.H

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
// 【方案 A 新增】isoAdvector 头文件
#include "isoAdvection.H"
// 【方案 B 新增】RDF 曲率 helper
#include "sigmaKRDF.H"
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
        "MPPICInterIsoFoam 衍生求解器：增加 runtime 可切换的 RDF 曲率选项。"
        "CSF 模式下退化到方案 A 行为；RDF 模式下窄带内用 RDF、窄带外 CSF。"
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

            // 【方案 B 新增】RDF 模式下更新 RDF 场
            //   share 关键点：此时 advector.surf() 已含最新 PLIC 重构
            //   （alphaEqnSubCycle 内部 advect() 刚刚重构过）
            //   CSF 模式下 rdfCurvPtr 为空，此分支 no-op。
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
