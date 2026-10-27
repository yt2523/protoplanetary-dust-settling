# 文献审查驱动的第一轮 Python 实验计划
日期：2026-10-27。状态：预先定义问题与验收；实验尚未实施/验证。
职责：Codex查文献、提出问题和判据；DeepSeek实现、测试和执行；Codex独立审核后提交重大改动。
基线：分支cpp-port/cpp17提交1b6fe48。先完成AUDIT_FIXES的验收门修复，再实施本计划。原仓库src/与原临时检出均不直接覆盖。

## 经核对的文献及适用范围
1. C. P. Dullemond & C. Dominik, 2004, A&A421,1075, DOI10.1051/0004-6361:20040284，arXiv astro-ph/0405226v1。Eq.21为固定、垂直等温气体中，沉降平流加气体密度加权的浓度扩散；无凝聚。本文提供结构依据，不保证本项目所有系数、盘密度和边界均相同。https://arxiv.org/html/astro-ph/0405226v1
2. A. N. Youdin & Y. Lithwick, 2007, Icarus192,588, DOI10.1016/j.icarus.2007.07.012，arXiv0707.2975v1。IV.2 Eq.18保留垂直惯性；Eq.21-28给垂向层厚关系。符号表的Sc_(x)为径向，tau_s=Omega*t_stop，St在该文另用t_stop/t_eddy，不能混用。本文不为直接在当前垂直方程添加1/(1+St²)提供充分依据。https://arxiv.org/html/0707.2975v1
3. R. Li & A. N. Youdin, 2021, ApJ，DOI10.3847/1538-4357/ac0e9f，arXiv2105.06042v2。强团聚门槛随颗粒大小、柱密度金属丰度、径向压力梯度和湍流水平变化；当前无反馈的一维沉降代码不能预测流致不稳定性或引力坍缩。https://arxiv.org/abs/2105.06042v2
4. 次要核查：Dullemond2004 Eq.4的停止时间前因子与当前sqrt(pi/8)写法不同，不据此认定项目错误；热速度/等温声速约定和Epstein近似前因子需单独交叉核验。本轮保留项目常数与停止时间。
HTML中的动态DRAFT日期不替代上述发表年与固定arXiv版本。未导入Zotero；引用与访问版本先在本项目保存。

## 问题与改变边界
P1：原dt=min(adv,diff)在组合显式算子上可能负密度，已有alpha=.1,d_lim=.7压力证据。
P2：原顶边界ghost浓度=0、gas=0，导致F_top,diff=alpha*rho_last/(2h)。常数浓度仅使内部扩散通量为0，不代表整个开放域的稳态。
P3：原层流z_real未跟随4H截断。先保留为基线限制，避免与P1/P2混合修正。
P4：any(rho_d/rho_g>1)仅是浓集标记，不能称成行星时间；新研究输出使用concentration_threshold，旧API保留兼容。
P5：声速裁剪不是终端速度近似有效性的证明。记录tau_s=Omega*t_stop、原始Mach与裁剪比例；tau_s≪1才是相关的小停止时间条件。本轮不改动阻力定律/Sc。

## 最小实现范围（独立Python模块，不替换C++默认基线）
建议research/settling_experiments.py与research/test_experiments.py；依赖仅已有NumPy/Matplotlib。
可继承当前原Python湍流类，继续调用其upwind_method。只提供两个显式实验开关：
- timestep=legacy或positive（正性联合步长）。
- top=vacuum_ghost（忠实原边界）或closed_total_flux（显式零总通量，属于边界模型改变）。
固定网格、气体、速度和加权扩散公式不变。closed可在调用原update后仅补回原顶通量在末格造成的损失，必须由通量测试验证。不要偷偷修层流网格或改Sc。
额外测试可直接调用隔离的扩散算子，令v=0；这只是检验分量，不是声称真实盘中没有重力。
比较同一网格、初值及实际无量纲结束时间。positive可用ceil(T/dt_bound)再全局dt=T/Nt精确对齐；标注这与legacy floor终止的区别，基线保持不变。不得把不同实际终止时间下的变化误称为算法误差。

## Codex独立代数推导（不是从文献直接抄出的算法）
g_i>0,q_i=rho_i/g_i，内部g_face=(g_left+g_right)/2。
内部对角损失率
lambda_i=[max(v_(i+1),0)+max(-v_i,0)]/h
         +alpha*(g_left_face+g_right_face)/(h²*g_i)。
下边界只保留右界面项。
顶边界vacuum的平流损失系数为abs(v_last)/h，扩散为alpha*(g_left_face+g_last/2)/(h²*g_last)；
closed顶边界去掉顶面通量项，保留朝向左界面的流出及扩散项。
对这个固定线性离散算子，非负初值保持非负的充分条件为dt*max(lambda)<=1；
用safety=.9留余量。不是新的物理规律或创新算法。
当alpha常数、v=-S*z/g、无速度裁剪、固定高斯g、零总通量时，直接积分F=v*rho-alpha*g*q'=0得到
rho_exact=C*g*exp[-(S/alpha)*(exp(z²/2)-1)]。
该解析平衡不能用于原vacuum顶边界，也不能用于被clip后的速度。
负速度原离散零通量平衡可递推
q_j=q_(j-1)/[1-h*v_j*g_j/(alpha*g_face_j)]，
再按h*sum(g*q)归一化。这是测试用离散平衡，不是替代时间演化的新算法。

## 预先固定的验证指标（执行前不得据结果放宽）
通用：每案例计划<=100000步，参考与实验各<=30秒。无参数大扫描，固定如下少量算例。
1. 正性对角系数检查：小网格N=8-16，用原update对基向量得到更新矩阵；对vacuum与closed分别核对推导的对角率及非负邻接系数，abs+rel门1e-12+1e-10*|reference|。
2. 已知压力案例：使用基线turb_diff_stress完全相同网格/物理输入和匹配实际结束时间。legacy必须复现显著负值；positive需全程有限、min>=-1e-14*initial_dust_peak，无clip修补，dt*lambda_max<=.9*(1+1e-12)。两者完整记录实际时间、质量、出流积分、各自min。质量残差<=1e-10*M_initial。旧负结果不是物理对照真值。
3. 稳定时间加密：positive用safety=.9/.45/.225在同一实际终止时间；前两者分别与最细结果比较，后一个RMS应小于前一个；报告RMS/Linf，不把数值时间加密误差当物理不确定度。
4. 扩散隔离、v=0、q=.01：closed演化100个安全步后仍满足初值逐点abs+rel门，且质量守恒；vacuum单步质量损失须匹配dt*alpha*rho_last/(2h)，并且非零。weighted不能换成对rho直接做普通Laplace。
5. 连续/离散平衡：N=32,64,128（最多256，不是大扫描）、L=2，S=.001,alpha=.01，确保未clip。用上述离散零通量递推剖面比较连续解析平衡，两者按同样离散质量归一化。RMS与Linf应随加密严格下降；记录观测阶（预期继承一阶迎风，不强行宣称精确阶数）。最细归一化RMS/连续峰值<.02。对离散平衡再实际演化若干安全步，变化满足abs+rel门，质量守恒。
6. 物理适用性诊断：记录tau_s范围、未clip/clip后Mach、cap_count、阈值事件；明确未模拟反作用/径向压力驱动/自引力。不以tau_s>.1之类任意阈值自动宣布模型无效，先作为启发式标记。额外全惯性ODE对照可留下一轮，不能声称本轮已做。

## 输出与Git关卡
research/RESULTS.md用中文列问题、文献位置、自己的推导、算例、判据、结果、未解决项；
research/results/ 保存CSV、JSON、图（稳态与压力对照，清楚区分时间/边界），以及命令/退出码、源码和参数哈希。
RUN_LOG追加每次重要动作及失败/修正。重大功能在Codex审核后单独提交并推送，不合并main，不改过去检查点。若测试未通过，保存失败证据和下一步，不写"已验证"。
## 物理独立审核补充（实验前，2026-10-27）
- N=8矩阵测试额外用dt=1.01/max(lambda)，确认至少一个对角系数为负；这不要求任意初值立刻负密度。
- vacuum与常数初始浓度在顶面不相容，首步通量随1/h放大，不能把这种跨网格变化直接叫实现错误。
- 边界敏感性使用alpha=.001、R100、N64、L2，在相同实际结束时间和相同安全时间步下比较vacuum/closed，隔离边界改变；报告质量与剖面差，不据阈值推出行星形成。
- 解析平衡仅适用于未clip负速度和零总通量。停止时间明确记tau_s=Omega*t_stop，避免文献St定义混淆。
- 方法验证不解决任意大粒径、反作用、真正SI或真实盘结构准确性；未验证项继续记录。