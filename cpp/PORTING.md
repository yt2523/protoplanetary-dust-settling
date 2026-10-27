# 移植映射与保留行为（PORTING）

## 来源与授权

- 源仓库：https://github.com/yt2523/protoplanetary-dust-settling
- 提交：`main` @ `6dd899a89203fb88c8724e69e10c208e4faf6d66`（只读检出，未做任何改动）
- 只读检出路径：`D:/Codex/tmp/phd-github-audit-2026-10-26/protoplanetary-dust-settling`
- **根目录未观察到 LICENSE 文件**，源码亦无许可头。因此本项目不臆造 MIT 或任何许可/授权声明，
  仅保留来源与作者署名：作者 `yt2523`，Imperial College London Year-3 Computational Physics (Project 2)。

## 类 → C++ 映射

| 原类 / 方法 | C++ 符号 |
| --- | --- |
| `PureAdvection1D_v2.__init__` | `dust::Advection` 构造 |
| `initial_function_points`（Gauss） | `Advection::initial` |
| `upwind_method`（有符号迎风） | `Advection::upwind` |
| `run_simulation` | `Advection::run` |
| `RealDiskDustSolver.__init__/calculate_natural_unit/choose_time_step/initial_dust_points/upwind_method/run_simulation/get_t_formation_years` | `dust::RealDiskDustSolver` |
| `DustSolver_Turbulence`（同上 + `alpha`、`v_cap_factor`、D=0 处理） | `dust::TurbulenceSolver` |
| 常量 `G,M_SUN,AU,Kb,M_H` | `dust_settling.cpp` 匿名命名空间常量 |

## 数值对齐的关键点（验收通过，严格门 EPS_ABS+EPS_REL·|P|）

- Euler 显式时间推进，`Nt = int(T/dt)`（向零截断），**无小数末步**；
  `t_actual = actual_steps * dt`（提前停止时 actual_steps < Nt，不再误用 `Nt*dt`）。
- 网格 `x/z = step*i`，`step = L/N`，等价 `np.linspace(..., endpoint=False)`。
- Gauss 初值 `sigma=0.07`、中心 `0.5`、峰值归一化到 `max`（**初始采样归一化**）；
  收敛用的 `gauss_exact` 现在也采用**初始采样归一化**（对未平移的采样峰值归一），以得到严格连续平移参考；
  原测试辅助函数对平移后峰值重新归一的差异已在 verify.py 注释中说明（偶数网格收敛不受影响）。
- 气体 `rho_g = exp(-0.5*z^2)`，尘埃初值 `rho_d = 0.01*rho_g`。
- `z` 域长度 `L = min(zmax*AU/H, 4.0)`（自然单位，上限 4）。
- 内界面迎风取 `v[i+1]`、气体密度取算术平均、梯度取 `rho_d/rho_g` 的差分（向量化加权扩散）。
- 下边界通量 = 0；上边界仅平流流出（`v<0` 时 `F_adv=0`）；上边界扩散 ghost `ratio=0`、`gas=0`
  （产生正的向上质量流失，未人为修正）。
- 成层判据为更新后 **`any(ratio>1)`**，而非中平面 3 格平均；`t_formation=(n+1)*dt`，`stop_when_cross` 提前 break。
- `v_cap_factor` 默认 **1**（裁剪到 c_s），湍流版专用；层流版无速度上限（见下基线怪癖）。

## 有意保留的基线怪癖（不做静默修正）

1. **层流 `z_real` 不裁剪**：`RealDiskDustSolver` 中 `z_real = linspace(0, real_z_max, Nz)`，
   而 `z`（自然单位）仍被 `L<=4` 裁剪；两者在 `zmax` 大到使 `L` 被钳制时不一致
   （`z_real_max ≠ real_z_max`，层流网格与 `z*H` 不一致是**唯一**的网格不一致点）。
2. **层流速度不设上限**：`RealDiskDustSolver` 不调用 `clip`，`v_real` 可任意大（湍流版才有 `v_cap_factor`）。
3. **湍流 `z_real = z*H`**：与层流的 `linspace` 不同，钳制后 `z_real_max = 4H < real_z_max`（L 上限失配）。
4. **D=0 时扩散 dt 为 +inf**：Python 中 `d_lim*h^2/(alpha*density_max)` 除零产生 `RuntimeWarning`
   并得到 `inf`，`dt=min(adv, inf)=adv`。C++ 显式写为 `max_dt_diff=+inf` 并输出等价告警。
5. `min(adv, diff)` 双限步长（`d_lim=0.7` 默认）并不保证合成算子的正性；按 PLAN 保留，
   主对比用 `adv=.3 diff=.2`，并保留一个小的有界压力用例（`adv=.3 diff=.7`，`alpha=0.1`）以
   **暴露而非掩盖**继承的负密度（`min_rho=-0.0124`），按 EXPECTED_INSTABILITY 单独报告，不裁剪。
6. 纯对流在 `abs(a)<1e-13` 时直接返回副本；`a` 的符号决定迎风方向；打印 "advective number … acceptable" 保留。

## 排除项（及原因）

- 教学类 `DustSolver1D`：`v(x)=-x` 的教学实现，非物理盘求解路径，PLAN 明确排除。
- 遗留标量 `legacy_turbulence_scalar.py`：仅为记录向量化优化历史而保留，非当前求解器。
- `PureAdvection1D`（方波父类）：被 `_v2`（Gauss）覆盖，不在验收范围。

## 已知限制（不声称）

- 仅对标量/教学数值身份做验收，不涉及物理结论或加速结论。
- C++ 二进制为 Linux ELF，经 `wsl.exe --exec` 以 argv 调用；Windows 侧无本机 C++ 编译器/CMake，
  CMake 4.4.4 由 `.tools/cmake` 提供并在 WSL 下实际运行（`build-cmake/`）。
