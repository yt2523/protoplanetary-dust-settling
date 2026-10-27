# 文献审查记录
访问日期：2026-10-27。Codex检索/审查，DeepSeek按research/PLAN.md实现测试。本文区分已检查正文与仅摘要筛选；不声称完成系统综述，也不把不同模拟假设的拟合式直接写入本模型。
## 已检查正文的模型基础
- Dullemond & Dominik (2004), A&A421,1075-1086，DOI10.1051/0004-6361:20040284；固定版本 https://arxiv.org/html/astro-ph/0405226v1 。检查2.1-2.2及3节，Eq.8/11/12/21。支持固定等温气体中的终端沉降与尘气比加权扩散结构；本项目盘密度、边界和停止时间前因子不因此自动获得验证。
- Youdin & Lithwick (2007), Icarus192,588-604，DOI10.1016/j.icarus.2007.07.012；固定版本 https://arxiv.org/html/0707.2975v1 。检查IV.2、附录A及符号表。垂向惯性、涡旋时间、层厚与径向Schmidt数分别讨论；该文tau_s=Omega*t_stop，而St=t_stop/t_eddy。速度裁剪不证明终端近似有效。
## 已核对摘要/元数据，尚未全面审查正文
- Li & Youdin (2021), Thresholds for Particle Clumping by the Streaming Instability，DOI10.3847/1538-4357/ac0e9f，https://arxiv.org/abs/2105.06042v2 。强团聚门槛依赖颗粒、柱密度金属丰度、压力梯度与湍流；当前any-cell>1标记不应称真实行星形成。
- Lim et al. (2024), ApJ969,130，DOI10.3847/1538-4357/ad47a2。读到的预印本是 https://arxiv.org/abs/2312.12508v1 （2023-12）；另由官方NSF论文存档核对2024发表信息。其3D、粒子自引力和外部强迫湍流结果不能无条件外推到本1D被动尘埃模型。未导入其拟合公式。
- Lim et al.，The Streaming Instability in 3D: Conditions for Strong Clumping，读版本 https://arxiv.org/abs/2509.18270v1 （2025-09-22）。摘要表明2D与3D门槛差异依停止时间变化；本条按预印本版本记录，不用网页之外的卷页替代实际读版本。
- Eriksson et al. (2026), Conditions for planetesimal formation via the streaming instability persist under turbulence driven by magnetorotational instability，https://arxiv.org/abs/2603.17195v2 （2026-08-19，摘要页标accepted A&A）。自洽MRI湍流与各向同性外强迫湍流的团聚影响不同；因此不能仅用一个扩散alpha推出所有湍流都阻止形成星子。本文仅筛选摘要，未复现MHD。
## 本轮决定与待办
1. 保留C++忠实基线，不改气体密度加权扩散，不直接加入径向Schmidt修正。
2. 先隔离检验组合显式步长的正性、原顶边界质量泄漏、零总通量平衡和实际终点对齐；这是数值方法验证，不是新物理发现。
3. 浓集阈值仅作诊断。无反作用、径向驱动、SI模态和自引力的当前代码不能单独确认星子形成。
4. 待核验的离散容量界见research/GRID_CAPACITY_REVIEW.md；不能把受网格容量限制的“未越阈值”解释成物理反证。
5. 原层流4H截断时物理网格不一致仍在基线中保留；未来单独修正版比较。停止时间前因子与声速/热速度约定需进一步交叉核查，未静默修改。
6. 原历史大扫描/观测约束未重跑；本次小算例不能确认原README的普适物理结论。
references.bib为可导入的书目元数据；本会话未进行Zotero写入。HTML中的动态DRAFT日期不替代论文发表年及固定arXiv版本。