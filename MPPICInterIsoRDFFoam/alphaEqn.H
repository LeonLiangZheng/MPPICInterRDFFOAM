// ============================================================================
// alphaEqn.H （方案 A 核心 —— 诊断变体 v3.1-noDAlphac）
// ----------------------------------------------------------------------------
// 目的：用 isoAdvector 几何输运替换 MULES，同时保持 MPPIC 的 αc 加权物理。
//
// 原 MPPIC α 方程（αc 加权）：
//     ∂(αc·α1)/∂t + ∇·(αc·α1·U) = Su + Sp·α1
//
// 完整方案 A 变形（链式法则 + 除以 αc）：
//     ∂α1/∂t + ∇·(α1·U) = Su/αc − (α1/αc)·[∂αc/∂t + ∇·(αc·U)]
//                       = Su/αc − (α1/αc)·dAlphac
//
// 即 Sp'·α1 + Su' 形式：
//     Su' = Su/αc
//     Sp' = Sp/αc − dAlphac/αc
//
// ─── 诊断变体修改（2026-05-19） ─────────────────────────────────────────
//   去掉 dAlphac/αc 注入项，目的：验证 RDF 内部 mottle 是否由 dAlphac
//   离散残差驱动 α1 偏离 1.0、触发 PLIC 在液滴内部产生假界面所致。
//
//   理论 dAlphac = ∂αc/∂t + ∇·(αc·U) ≡ 0（载流体连续性方程），
//   但本求解器中：
//     - αc 由 kinematicCloud.theta() 粒子统计投影得到（非守恒 PDE 解）
//     - phi 由动量方程解出的 U 插值得到
//   两者不被同一守恒律耦合，dAlphac 在粒子运动 cell 处出现量级
//   ~10⁵–10⁶ s⁻¹ 的离散残差 spike。注入 α 方程后驱使 α1 偏离 1.0
//   一个 ~0.1 量级的扰动，远超 surfCellTol（默认 1e-6），触发 PLIC
//   在液滴内部把粒子 cell 误判为 partial cell → 假界面 → 污染 RDF。
//
//   代价：牺牲严格的 αc·α 离散守恒等价性。
//   预期影响：稀粒子云 (αc > 0.95) 下水质量漂移 << 1%，可接受。
//   日后若需恢复严格守恒，考虑：dAlphac 加界面 mask、或粒子 cell 限幅。
// ────────────────────────────────────────────────────────────────────────
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

    // --- 步骤 3：变形 Sp、Su（诊断变体：去掉 dAlphac 注入项）
    //    完整方案 A 应该是：
    //        Su' = Su/αc
    //        Sp' = Sp/αc − dAlphac/αc
    //    现在仅保留 Su/αc 与 Sp/αc，dAlphac 项被有意删除。
    //    详见文件头部注释段。
    Su = Su/alphacSafe.internalField();
    Sp = Sp/alphacSafe.internalField();

    // --- 步骤 4：isoAdvector 几何输运
    //    内部会调用 surf().reconstruct() 做 PLIC，然后按 dt 解析几何
    //    通量。α1 被原地更新，advector.alphaPhi() 存储 α1 的面通量
    //    (α1f·phi，量纲 m^3/s，不含 αc 加权)。
    advector.advect(Sp, Su);

    // 网格运动修正还原
    if (mesh.moving())
    {
        U += fvc::reconstruct(mesh.phi());
    }

    // --- 步骤 5：组装 αc 加权的质量通量 rhoPhi
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
