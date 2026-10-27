# 原行星盘尘埃沉降求解器 —— C++17 忠实移植

本项目是 [protoplanetary-dust-settling](https://github.com/yt2523/protoplanetary-dust-settling)
（`main` @ `6dd899a89203fb88c8724e69e10c208e4faf6d66`）三个核心求解器的 **C++17 忠实移植**，
仅用标准库实现，不依赖 NumPy。移植目的为数值行为对照的验收验证，不提出任何加速或物理新发现。

## 范围

| 原 Python 模块 / 类 | 移植状态 |
| --- | --- |
| `pure_advection.py` / `PureAdvection1D_v2`（Gauss 初值、有符号、周期、一阶迎风） | ✅ |
| `dust_solver.py` / `RealDiskDustSolver`（层流、标量循环、速度**不裁剪**） | ✅ |
| `dust_solver_turbulence.py` / `DustSolver_Turbulence`（向量化加权扩散、速度**裁剪**） | ✅ |
| 教学用 `DustSolver1D` 与遗留标量 `legacy_turbulence_scalar.py` | ❌ 有意排除（见 PORTING.md） |

## 构建

C++ 二进制由 D 盘后端 Ubuntu WSL 的 g++ 15.2.0 + **CMake 4.4.4**（`.tools/cmake` 内自带，非系统安装）编译。

**标准命令（读者可在装有 CMake 的任意环境使用）：**

```bash
cmake -S . -B build-cmake -DCMAKE_BUILD_TYPE=Release
cmake --build build-cmake -j4
```

**本仓库实际使用的本地自带 CMake 绝对命令（Windows WSL 后端）：**

```bash
wsl.exe bash -lc 'cd /mnt/d/Codex/tmp/dust-cpp-publish-20261027/cpp && \
  TMPDIR=/mnt/d/Codex/tmp \
  /mnt/d/Codex/downloads/protoplanetary-dust-settling-cpp/.tools/cmake/cmake/data/bin/cmake \
    -S . -B build-cmake -DCMAKE_BUILD_TYPE=Release && \
  /mnt/d/Codex/downloads/protoplanetary-dust-settling-cpp/.tools/cmake/cmake/data/bin/cmake \
    --build build-cmake -j4'
```

> 产物路径 `build-cmake/dust_settling`（本出版布局内）；路径由 `verify.py`/`regression.py`
> 从自身 `ROOT` 或 `--binary` 派生并做 WSL `/mnt/<drive>/...` 转换，无硬编码绝对交付路径。

`CMakeLists.txt` 配置静态库 `dust_settling_lib` + CLI `dust_settling`，并显式禁用浮点收缩
（GNU/Clang 用 `-ffp-contract=off`，MSVC 用 `/fp:strict`），无 fast-math。
Windows 侧无本机 C++ 编译器/CMake，故全部编译在 WSL 下完成，产物为 Linux ELF。

## CLI 用法

```text
dust_settling advection   v dt T Nx outdir
dust_settling laminar     total_time Nz a_lim R_AU a_micron zmax dust_density stop outdir
dust_settling turbulence  total_time Nz a_lim d_lim R_AU a_micron zmax D dust_density v_cap stop outdir
```

CLI 对参数做严格校验：各模式精确 `argc`、doubles 非 NaN/inf、正整数（N 界 `[2,1000000]`）、
物理量正负号约束、步数在 `int` 转换前先做有限/范围检查（`checked_step_count`，在 run 循环前阻断，
拒绝 `1e308/1e-308` 这类巨量步、不含运行后诊断兜底）、异常返回非零、写文件失败返回非零
（开启 stream 异常并检查 open/final-flush）。每个运行在 `outdir` 写入：
- `profile.csv` —— 严格 4 个逗号字段，列头不含逗号/单位歧义
  （`z_over_H,rho_d_over_rho_g0,rho_g_over_rho_g0,z_coordinate_m`；
    advection 为 `x_over_L,rho_final_normalized,rho_initial_normalized,x_over_L_aux`）
- `diagnostics.txt` —— `key value` 标量诊断（dt、Nt、actual_steps、t_actual、t_actual_years、
  质量、边界通量、min/max、形成时间、cap_count 等）

## 验收

```bash
C:/Users/PC/AppData/Local/Programs/Python/Python312/python.exe -B verify/verify.py
```

该脚本 `import` 原始 Python 模块作为基准（`--source-repo` / `--binary` 可覆盖默认路径），
经 `wsl.exe --exec` 以 argv 直接调 C++ 二进制，逐点比较。每个 Python 参考用例在**受限子进程**
（30 s 超时、无重试，`verify/_ref_worker.py`）中隔离运行。结果写入 `validation/results.csv` 与
`validation/results.json`，并为每个用例保存 Python 参考 profile。本轮（见 `验证报告.md`）：
**17 项数值校验 PASS + 1 项 EXPECTED_INSTABILITY（继承的负密度，单独报告不裁剪）+ 1 项
convergence PASS + 1 项 branch-coverage PASS**，0 失败；密度/坐标逐点差 ≤ EPS_ABS + EPS_REL·|P|，
质量平衡残差按独立积分质量 + 出流积分通量复核（不信任发出的 mass_balance_error），收敛阶 0.9003。
一个 diff=.7 压力用例暴露继承的负密度，按 EXPECTED_INSTABILITY 单独报告。

**小型运行示例（单个用例 + 绘图，均使用已保存 profile，不重跑求解器）：**

```bash
# 单次 advection 运行（写 outdir/profile.csv + diagnostics.txt）
wsl.exe --exec build-cmake/dust_settling advection 1 0.01 0.5 50 /mnt/d/<path>/outdir

# 轻量绘图（读取已保存的 C++/Python profile，不重跑求解器）
C:/Users/PC/AppData/Local/Programs/Python/Python312/python.exe -B verify/plot.py
```

**CLI 稳健性回归（20 项，含负 v_cap 拒绝、来源哈希前后一致性）：**

```bash
C:/Users/PC/AppData/Local/Programs/Python/Python312/python.exe -B verify/regression.py
```

**验收门单元测试（见证质量/最小值变异从不被误判为 PASS）：**

```bash
C:/Users/PC/AppData/Local/Programs/Python/Python312/python.exe -m pytest verify/test_gates.py -v
```

来源 SHA 前后一致、精确命令与结果见 `validation/provenance.json`。

## 目录

```
src/              C++17 库（dust_settling.hpp / .cpp）与 CLI（main.cpp）
verify/verify.py  验收脚本（TDD：先验证缺失二进制失败，再验证数值一致）
validation/       机器可读证据（CSV/JSON、各用例诊断）与中文报告
CMakeLists.txt    CMake 构建配置（库 + CLI）
PORTING.md        移植映射、保留的怪癖、限制与来源声明
```

## 许可与署名

原始仓库根目录**未观察到 LICENSE 文件**，源文件亦无许可头，因此本项目**不臆造 MIT 或任何授权**，
仅保留来源与作者署名：原仓库 `yt2523/protoplanetary-dust-settling`，提交 `6dd899a…`，
Imperial College London Year-3 Computational Physics (Project 2)。详见 PORTING.md。
