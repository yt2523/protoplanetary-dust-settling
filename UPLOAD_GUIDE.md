> 这份文件是给你自己看的说明，上传前删掉（留着也不影响别人读代码）。

# 这个仓库是怎么整理出来的

## 你的原始文件动过吗？

**没有。** `C:\Users\PC\OneDrive - ...\work\CP_project` 和 `D:\学学学学\作业\CP project`  
两个目录都是只读访问，一个字节没改。这里是一份重新组织过的副本。

## 和 c1 项目不一样：这次没有大文件问题

c1 那个项目最大的坑是 1.5 GB 的 CSV 顶穿 GitHub 的 100 MB 单文件上限。  
这个项目**完全没这个问题** —— 代码目录只有 225 KB，因为它是纯计算项目，  
所有数据都是运行时算出来的，没有需要存盘的观测数据集。整个仓库加上图和报告 PDF 也才 1 MB 出头。

所以这次的重点不是"删什么"，而是"怎么让人看懂"。

## 具体做了五件事

### 1. 分开 src/ 和 tests/

原来 10 个 `.py` 平铺在一个目录里，求解器和跑图脚本混在一起，看的人分不清主次。现在：

- `src/` —— 四个求解器模块，这才是你要展示的工程能力
- `tests/` —— 五个驱动脚本，每个对应报告里的一节

因为分了目录，测试脚本 import 不到求解器了，所以加了 `tests/_path_setup.py`  
把 `../src` 注进 `sys.path`，五个测试文件开头各加一行 `import _path_setup`。  
**求解器本身的算法逻辑一行没改。**

### 2. 文件名从"第几节"改成"做什么"

这是最重要的一步。原来的名字只有你自己懂：

| 原名                              | 新名                              | 对应报告     |
| ------------------------------- | ------------------------------- | -------- |
| `test_pure_advection.py`        | `test_convergence_advection.py` | Fig 1, 2 |
| `test_r.py`                     | `test_numerical_stability.py`   | Table 1  |
| `test_real_problem_section4.py` | `test_settling_vs_radius.py`    | Fig 3    |
| `test_section4.1.py`            | `test_turbulence_alpha_scan.py` | Fig 4–6  |
| `test_section5.py`              | `test_grain_size_scan.py`       | 尺寸扫描     |

`test_section4.1.py` 对招聘方是零信息量 —— 他不知道你报告第 4.1 节写了什么。

### 3. 处理两个版本的湍流求解器

你有 `dust_solver_turbulence.py`（标量循环）和 `dust_solver_turbulence_vector.py`（向量化）。  
向量化版才是主力，但名字带 `_vector` 后缀反而像个附属实验品。所以：

- 向量化版 → `src/dust_solver_turbulence.py`（**主求解器**）
- 标量版 → `src/legacy_turbulence_scalar.py`（保留作对比）

**没删旧版** —— "先写朴素循环，发现太慢，重构成 NumPy 切片"这个过程本身是加分项，  
README 专门写了一段。

### 4. 用你论文里的真图

我原本想用你的代码重新跑图，但 α=5×10⁻³ 那档的抛物型限制 Δt ∝ h²/α 把步长压得极小，  
5 Myr 积分要几十分钟。既然你论文里的图就是同一套代码跑出来的，直接用更准也更快：

| 论文原名           | 仓库里的名字                                   |
| -------------- | ---------------------------------------- |
| `Figure_1.png` | `results/advection_profiles.png`         |
| `Figure_2.png` | `results/advection_convergence.png`      |
| `Figure_7.png` | `results/settling_time_vs_radius.png`    |
| `Figure_4.png` | `results/turbulence_alpha_scan.png`      |
| `Figure_5.png` | `results/turbulence_alpha_scan_weak.png` |
| `Figure_6.png` | `results/turbulence_profiles_scaled.png` |

`Figure_3.png` 没用 —— 它的纵轴范围是 8777~8786 yr，和 `Figure_7` 的 1–2 Myr 差三个量级，  
应该是早期版本或不同参数下的输出。我用了 `Figure_7`，因为它带完整拟合和置信带，  
而且数值和你论文正文写的"1–2 Myr"一致。**如果 Figure_3 其实是你要的那张，告诉我，我换回来。**

报告 PDF 也放进 `docs/report.pdf` 了，别人想看细节可以直接翻。

### 5. 加了四个 GitHub 必需文件

- `README.md` —— 主文档，数字全部来自你的报告
- `results/RESULTS.md` —— 结果速查表，不想跑代码也能看到结论
- `.gitignore` —— 屏蔽 `__pycache__/`、`.vscode/`、`*.zip`、`*.docx`
- `requirements.txt` —— 只有 numpy 和 matplotlib

---

## ⚠️ 有三处我先前写错了，已按你论文改正

我在读到你的报告之前，凭代码推断写过一版 README，其中三处结论是错的。现在都改了，  
但你自己也要知道正确的说法是什么，面试被问到别说反：

**① 沉降时间和半径的关系**

- ❌ 我原先写："内盘密度高沉降快，往外逐渐变慢"
- ✅ 你论文的真实结论：**20 AU 约 1.8–2.0 Myr，100 AU 约 0.99 Myr，超过 50 AU 后半径依赖基本可忽略**

这个才是有意思的点：往外走气体变稀（沉降变快），但也变冷、引力束缚变弱（沉降变慢），  
两个效应几乎抵消。这是"看起来平淡的结果背后有物理"，比单调趋势更值得讲。

**② 微米颗粒能不能达到阈值**

- ❌ 我原先写："α ≳ 10⁻³ 时达不到阈值"（暗示小 α 时能达到）
- ✅ 真实结论：**所有 α 都达不到**。中平面尘气比只有 0.043~0.25，全部远低于 1

**③ 大颗粒能不能靠长大解决**

- ❌ 我原先写："厘米级以上颗粒能在盘寿命内穿越阈值"
- ✅ 真实结论：**1 Myr 内没有任何尺寸（直到 1 cm）能穿越**。因为盘太密，阻力太强，  
  沉降速度随尺寸增长得比稀盘慢得多

这三条其实让你的项目更有价值 —— 它是一个**有明确负面结论的定量研究**，  
而且你还解释清楚了为什么（dense disk 的阻力效应），并对比了稀盘该是什么行为。  
比"我复现了教科书结论"强。

---

## 上传步骤

```bash
cd "D:/workbuddy/github/protoplanetary-dust-settling"

git init
git add .
git commit -m "1-D advection-diffusion solver for dust settling in protoplanetary disks"

# 先在 GitHub 网页建空仓库，名字建议 protoplanetary-dust-settling
git remote add origin https://github.com/ty12312322/protoplanetary-dust-settling.git
git branch -M main
git push -u origin main
```

推之前确认没有超标文件（应该是空输出）：

```bash
find . -size +50M -not -path "./.git/*"
```



---

## 简历上怎么写这条

你现在简历里"PDE 数值求解器开发：长时域物理演化模拟"这个标题太抽象，看不出难度。建议换成：

> **原行星盘尘埃沉降 PDE 求解器** — 从零实现一维 advection–diffusion 显式有限差分求解器，  
> 用于百万年尺度尘埃垂直演化模拟。对流项用逐界面迎风保证正定性与单调性，扩散项按  
> ρ_g ∂z(ρ_d/ρ_g) 形式离散使良好混合态成为格式真 steady state；自适应时间步同时满足  
> 双曲 CFL 与抛物稳定条件，并对超声速区加软速度上限避免步长退化到 5×10⁻⁶ yr。  
> 以自然单位无量纲化消除量纲系数的数量级跨度，使 Myr 级积分可行。在解析解上验证一阶收敛，  
> 并通过网格 / CFL / 速度上限三组扫描将数值不确定度界定在 3%（t_settle 对 CFL 变化 <0.03%）。  
> 将标量循环重构为 NumPy 向量化，使多参数扫描可负担。定量结论：该致密盘模型下  
> 1 Myr 内无任何颗粒尺寸可达星子形成阈值，需同时满足颗粒生长与弱湍流（α≈10⁻⁵–5×10⁻⁵）。  
> <https://github.com/ty12312322/protoplanetary-dust-settling>

关键词密度是重点，招聘方会认这些：**有限差分、迎风格式、正定性/单调性、CFL 与抛物稳定性条件、  
无量纲化、收敛阶验证、数值不确定度量化、向量化性能优化**。

这些直接对口你的**第三梯队（CAE / EDA 数值计算）**，也是**第二梯队（仿真测试）**&#x770B;重的  
"懂求解器内部而不只是调 API"的硬证据。

短版（一行，用于篇幅紧的简历）：

> 一维 advection–diffusion 有限差分求解器（迎风+中心差分，双稳定性条件自适应步长，  
> 无量纲化，收敛阶验证，数值不确定度 3%），NumPy 向量化后支持 Myr 级多参数扫描。

配合 c1 那个项目，你 GitHub 上就同时有了「数值方法 / 物理建模」和  
「机器学习 / 统计推断」两条线，覆盖面比只有一个项目强得多。
