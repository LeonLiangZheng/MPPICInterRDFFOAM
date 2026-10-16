// ============================================================================
// alphaEqn.H （方案 A 核心）
// ----------------------------------------------------------------------------
// 目的：用 isoAdvector 几何输运替换 MULES，同时保持 MPPIC 的 αc 加权物理。
//
// 原 MPPIC α 方程（αc 加权）：
//     ∂(αc·α1)/∂t + ∇·(αc·α1·U) = Su + Sp·α1
//
// 方案 A 变形（链式法则 + 除以 αc）：
//     ∂α1/∂t + ∇·(α1·U) = Su/αc − (α1/αc)·[∂αc/∂t + ∇·(αc·U)]
//
// 等价写成 Sp'·α1 + Su' 形式：
//     Su' = Su/αc
//     Sp' = Sp/αc − dAlphac/αc ,    其中 dAlphac = ∂αc/∂t + ∇·(αc·U)
//
// 左端正是 isoAdvector::advect(Sp', Su') 可直接消化的形式。
// ============================================================================

{
    // 网格运动修正（对 MPPIC 多相场景通常 mesh 不动，但保留兼容）
    if (mesh.moving())
    {
        U -= fvc::reconstruct(mesh.phi());
    }

    // --- 步骤 1：先按标准 alphaSuSp.H 初始化 Sp、Su（吸收 fvOptions 源项）
    #include "alphaSuSp.H"

    // --- 步骤 2：αc clipping，防止极端接近零时后续除法爆炸
    //    【方案 A 变形点】只影响 Sp'/Su' 的分母，不改变 α1 输运
    const dimensionedScalar alphacMinClip
    (
        "alphacMinClip", dimless, alphacMinClipValue
    );
    const volScalarField alphacSafe("alphacSafe", max(alphac, alphacMinClip));

    // --- 步骤 3：计算 αc 的变化项 dAlphac = ∂αc/∂t + ∇·(αc·U)
    //    量纲：dimless/dimTime （与 Sp 匹配）
    //    注：用 alphacf*phi 代替 αc·U，保持面中心有限体积离散
    //    【方案 A 变形点】右端"αc 链式修正项"
    const volScalarField dAlphac
    (
        "dAlphac",
        fvc::ddt(alphac) + fvc::div(alphacf*phi)
    );

    // --- 步骤 4：变形 Sp、Su
    //    Sp' = Sp/αc − dAlphac/αc
    //    Su' = Su/αc
    //    【方案 A 变形点】核心一行（Internal 运算）
    Su = Su/alphacSafe.internalField();
    Sp = Sp/alphacSafe.internalField()
       - dAlphac.internalField()/alphacSafe.internalField();

    // --- 步骤 5：isoAdvector 几何输运
    //    内部会调用 surf().reconstruct() 做 PLIC，然后按 dt 解析几何
    //    通量。α1 被原地更新，advector.alphaPhi() 存储 α1 的面通量
    //    (α1f·phi，量纲 m^3/s，不含 αc 加权)。
    advector.advect(Sp, Su);

    // 网格运动修正还原
    if (mesh.moving())
    {
        U += fvc::reconstruct(mesh.phi());
    }

    // --- 步骤 6：组装 αc 加权的质量通量 rhoPhi
    //    【接力文档 §7 修正点 —— 重要】
    //    原文档写的 "rhoPhi = advector.alphaPhi()*(rho1f-rho2f) + alphacf*phi*rho2f"
    //    缺少第一项的 αcf 加权，导致与原 UEqn/pEqn 的 αc 加权物理不一致。
    //
    //    推导：
    //      rhoPhi 应等于 αcf·phi·ρ_mix_f
    //                 = αcf·phi·[α1f·ρ1 + (1−α1f)·ρ2]
    //                 = αcf·α1f·phi·(ρ1−ρ2) + αcf·phi·ρ2
    //                 = αcf · advector.alphaPhi() · (ρ1f−ρ2f) + αcf·phi·ρ2f
    //
    //    对照原 MPPICInterFoam：
    //      alphaPhi10_MULES ≡ αcf·α1f·phi  (phiCN = αcf·phi → fvc::flux 的结果)
    //      rhoPhi_MULES = alphaPhi10_MULES*(ρ1f−ρ2f) + phiCN*ρ2f
    //    完全对应。
    {
        #include "rhofs.H"

        // 同时维护 alphaPhi10 字段，使 alphaCourantNo.H 等外部代码继续可用
        alphaPhi10 = alphacf * advector.alphaPhi();

        rhoPhi = alphaPhi10*(rho1f - rho2f) + alphacf*phi*rho2f;
    }

    alpha2 = 1.0 - alpha1;
    mixture.correct();
}

Info<< "Phase-1 volume fraction = "
    << alpha1.weightedAverage(mesh.Vsc()).value()
    << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
    << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
    << endl;
