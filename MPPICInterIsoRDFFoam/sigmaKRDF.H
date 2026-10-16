/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::sigmaKRDF

Description
    基于 RDF (Reconstructed Distance Function) 的 σ·κ 曲率计算 helper。

    核心思路：
      - 窄带内：κ = -div(grad(RDF)/|grad(RDF)|)  —— RDF 对 parasitic current 更友好
      - 窄带外：κ = mixture.K()                  —— CSF，保护 Kistler 接触角 BC
      - 连续 mask 由 rdf.cellDistLevel() 提供线性渐变（不是 boolean 硬阶跃）
      - κ 层 blend，外乘 σ

    架构决策（share 模式）：
      - 不自建 isoAdvection，仅持主 advector 的 const 引用
      - update() 读取 advector.surf() 的 PLIC 重构结果（centre/normal/interfaceCell），
        避免重复重构（64 rank HPC 上这次开销非零）
      - 调用时序：alphaEqnSubCycle.H 之后 → rdfCurv.update() → UEqn.H

    v2 清理 vs 旧 sigmaKRDF 的 5 处差异：
      1) 直接 mixture.sigma()，不手工解析 transportProperties
      2) 默认 nSmoothIter = 0（保护 Kistler 接触角动力学）
      3) 连续 mask 用 cellDistLevel，不是 boolean + Laplacian 磨
      4) deltaN 跟网格缩放：1e-8/cbrt(meanCellVolume)
      5) κ 层 blend（而非 σκ 层），外乘 σ —— 数学等价但更自然

\*---------------------------------------------------------------------------*/

#ifndef sigmaKRDF_H
#define sigmaKRDF_H

#include "fvCFD.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "isoAdvection.H"
#include "reconstructedDistanceFunction.H"

namespace Foam
{

class sigmaKRDF
{
private:

    // Private data

        const fvMesh& mesh_;

        //- 持 const 引用：只读 sigma、K（CSF 曲率）
        const immiscibleIncompressibleTwoPhaseMixture& mixture_;

        //- share 模式：持主 advector const 引用，复用 surf() 的 PLIC 重构结果
        // const isoAdvection& advector_;
        isoAdvection& advector_;

        //- RDF 场（继承 volScalarField）
        autoPtr<reconstructedDistanceFunction> rdfPtr_;

        //- 窄带宽度（默认 4 层 cell）
        const label bandWidth_;

        //- 法向平滑迭代次数（默认 0，保护 Kistler 接触线尖锐曲率）
        const label nSmoothIter_;

        //- 梯度稳定项，dimless（grad(RDF) 本身 dimless，因 RDF 是 dimLength）
        const dimensionedScalar deltaN_;
        
        //- 表面张力 σ（从 PIMPLE 子字典的 RDFsigma 参数读，N/m）
        //  注：v2212 的 immiscibleIncompressibleTwoPhaseMixture 不公开 sigma()，
        //  这里避开 interfaceProperties 的访问限制，显式由 fvSolution 提供
        const dimensionedScalar sigma_;

    // 禁止拷贝
    sigmaKRDF(const sigmaKRDF&) = delete;
    void operator=(const sigmaKRDF&) = delete;


public:

    // Constructor
    sigmaKRDF
    (
        const fvMesh& mesh,
        const immiscibleIncompressibleTwoPhaseMixture& mixture,
        //const isoAdvection& advector,
        isoAdvection& advector,
        const dictionary& dict
    );


    // Public member functions

        //- 更新 RDF 场；前提是 advector.surf() 已含最新 PLIC 重构
        void update();

        //- 返回 σ·κ 场（窄带内 RDF + 窄带外 CSF 的 blend）
        tmp<volScalarField> sigmaK() const;

        //- 访问底层 RDF（调试/可视化）
        const reconstructedDistanceFunction& rdf() const { return rdfPtr_(); }
};

} // End namespace Foam

#endif

// ************************************************************************* //
