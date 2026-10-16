// ============================================================================
// alphaControls.H （方案 A 新版）
// ----------------------------------------------------------------------------
// 读取 fvSolution 的 solvers["alpha.water.*"] 子字典。
// 由于 MULES 整个被 isoAdvector 替换，原 MPPICInterFoam 的
//   MULESCorr / nAlphaCorr / alphaApplyPrevCorr / icAlpha / scAlpha
// 等控制项全部弃用。
// 只保留：
//   - nAlphaSubCycles  (isoAdvector subCycle 次数)
//   - alphacMinClip    (αc 的下限 clip 值，防 αc→0 时除法爆炸)
// ============================================================================

const dictionary& alphaControls = mesh.solverDict(alpha1.name());

label nAlphaSubCycles
(
    alphaControls.get<label>("nAlphaSubCycles")
);

// αc clipping 参数：
//   - Leon 实测 Min(alphac) ≈ 0.5
//   - 保守估计未输出时步可能跌至 ≈ 0.1
//   - 默认 0.05：远低于实际工况、仅作为 αc 异常接近零时的防 NaN 安全阀
// case 里可通过 alphacMinClip 覆盖默认值
const scalar alphacMinClipValue
(
    alphaControls.getOrDefault<scalar>("alphacMinClip", 0.05)
);
