# 文献驱动的 Python 数值实验
本目录与 C++ 忠实移植分开。DeepSeek编写实验和测试，Codex审查文献、推导、结果与提交；原仓库src和C++默认方程未修改。包装器调用原湍流类的upwind_method，并非另一套普通密度扩散方程。

从完整分支的cpp/research运行（需要已有NumPy、Matplotlib、pytest）：

    python -B run_final_checks.py

独立D盘下载目录先在PowerShell指定原Python的src目录：

    $env:TEMP='D:\Codex\tmp'
    $env:TMP=$env:TEMP
    $env:PYTHONDONTWRITEBYTECODE='1'
    $env:DUST_SOLVER_SRC='D:\Codex\tmp\phd-github-audit-2026-10-26\protoplanetary-dust-settling\src'
    python -B run_final_checks.py

源路径按显式参数/环境变量、完整分支父级src、本机原临时克隆选择；结果记录实际路径及文件SHA256。

- settling_experiments.py：TurbulenceExperiment，两个独立开关。timestep=legacy保留原步长及floor终点；positive用联合算子正性步长，缺省在请求T结束，比较时显式传target_time_natural=legacy_actual_time。top=closed_total_flux改变物理边界；vacuum_ghost保留原边界。无负密度裁剪。
- test_experiments.py：35项矩阵、扩散、质量、平衡、来源及时间对齐检查。
- test_finalization.py：2项真实诊断函数和已保存绘图数据回归。
- run_experiments.py：固定少量算例、CSV/JSON和PNG；图从保存剖面读取。密度为rho_d/rho_g0，坐标为z/H，时间含无量纲值与年。
- audit_grid_capacity.py：仅原层流初始化和已有剖面读数，无时间积分；核验尾部质量约束和网格容量界。不是通用物理API。
- run_final_checks.py：依次执行实验、37项测试、容量审计，每子进程30秒超时，保存实际返回码和stdout/stderr。
- PLAN.md：运行前判据；RESULTS.md：已审核结果；RUN_LOG.md：失败、修正与执行分工；GRID_CAPACITY_REVIEW.md：适用条件及推导。
- results/：当前机器可读证据和图。以final_check_summary.json、final_*_stdout.txt为本次执行证据。

文献阅读深度与书目见上级LITERATURE_AUDIT.md、references.bib；没有进行Zotero写入。项目不是含反馈、自引力或SI的星子形成模拟。

限制：固定气体/速度、常alpha；连续平衡只适用未裁剪速度和零总通量。每次运行新建对象；零算子、强制步数/时间步及任意参数组合未作通用API验收。循环内时间预算每500步检查，外层30秒超时才是本轮的硬限制。未处理原层流4H截断后物理坐标不一致；未验证全惯性阻力模型或任意粒径物理适用性。