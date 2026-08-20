import _path_setup  # noqa: F401  (adds ../src to sys.path)
import numpy as np
import matplotlib.pyplot as plt
from pure_advection import PureAdvection1D_v2 as PureAdvection1D

# analytic Gaussian, periodic boundary
def gauss_exact(x, v, t, L=1.0, x0=0.5, sigma=0.07):
    # back-traced position with periodic wrap
    x0_back = (x - v * t) % L
    y = np.exp(-0.5 * ((x0_back - x0) / sigma)**2)
    # normalise peak to 1
    y /= y.max()
    return y

def err(y_num, y_ex):
    d   = y_num - y_ex
    l2  = np.sqrt(np.mean(d**2))
    lin = np.max(np.abs(d))
    return l2, lin


def all_plot():
    import numpy as np
    import matplotlib.pyplot as plt

    # ----- 字号整体调大 -----
    plt.rcParams.update({
        "font.size": 14,             # 主字号
        "axes.titlesize": 16,        # 标题字号
        "axes.labelsize": 14,        # x, y 轴标签
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 12
    })

    v   = 1.0
    cfl = 0.5
    L   = 1.0
    T   = 0.5

    nx0 = 500
    h0  = L / nx0
    dt0 = cfl * h0 / abs(v)

    sim0 = PureAdvection1D(
        c_velocity=v,
        time_step=dt0,
        total_time=T,
        spatial_step=nx0
    )
    y0_num = sim0.run_simulation()
    x0     = sim0.x
    t_end  = sim0.Nt * sim0.dt
    y_ex0  = gauss_exact(x0, v, t_end)

    # initial (Gaussian at t=0)
    y_init = gauss_exact(x0, v, 0.0)

    plt.figure(figsize=(12,5), dpi=200)
    plt.step(x0, y_init, where="mid", alpha=0.5, label="initial")
    plt.step(x0, y_ex0, where="mid", linestyle="--", color="k", label="analytic")

    # Nx set
    nx_set = [50, 150, 200, 500]
    cols   = ["tab:purple", "tab:red", "tab:orange", "tab:green"]

    for nx, cc in zip(nx_set, cols):
        h  = L / nx
        dt = cfl * h / abs(v)

        sim = PureAdvection1D(
            c_velocity   = v,
            time_step    = dt,
            total_time   = T,
            spatial_step = nx
        )
        y_num = sim.run_simulation()
        x     = sim.x
        t_end = sim.Nt * sim.dt

        y_ex  = gauss_exact(x, v, t_end)
        rms, mmax = err(y_num, y_ex)
        lb = f"Nx={nx} (rms={rms:.3e}, max={mmax:.3e})"

        plt.step(x, y_num, where="mid", color=cc, alpha=0.85, label=lb)

    plt.xlabel("x")
    plt.ylabel("ρ(x)")
    plt.title(f"Pure advection of a Gaussian pulse at t={t_end:.3f}")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.show()


# ------- 分辨率对比图（原来的第一个 conv_plot，我只改了名字+高斯） -------
def resol_plot():
    v   = 1.0
    cfl = 0.5
    L   = 1.0
    T   = 0.5

    # pick 4 resolutions to show visually
    nx_set = [50, 150, 200, 500]
    colors = ["tab:purple", "tab:red", "tab:orange", "tab:green"]
    lbls   = ["Nx=50", "Nx=150", "Nx=200", "Nx=500"]

    plt.figure(figsize=(8,4))

    for nx, cc, lb in zip(nx_set, colors, lbls):
        h  = L / nx
        dt = cfl * h / abs(v)

        sim = PureAdvection1D(
            c_velocity   = v,
            time_step    = dt,
            total_time   = T,
            spatial_step = nx
        )
        y_num = sim.run_simulation()
        x     = sim.x
        t_end = sim.Nt * sim.dt

        y_ex = gauss_exact(x, v, t_end)

        plt.step(x, y_num, where="mid", color=cc, alpha=0.8, label=lb)

    # analytic curve
    y_ex_plot = gauss_exact(x, v, t_end)
    plt.step(x, y_ex_plot, where="mid", linestyle="--", color="k", label="analytic")

    plt.xlabel("x")
    plt.ylabel("Density ρ(x)")
    plt.title("Resolution comparison of upwind advection (Gaussian pulse)")
    plt.legend()
    plt.tight_layout()
    plt.show()


# ------- 真正的收敛 log–log 图（保持原结构，只换成高斯） -------
def conv_plot():
    import numpy as np
    import matplotlib.pyplot as plt

    v   = 1.0
    cfl = 0.5
    L   = 1.0
    T   = 0.5

    nx_list = list(range(50, 501, 50))   # 50,100,...,500
    h_list  = []
    e_list  = []

    print("\nConvergence test (Gaussian pulse, CFL=0.5):")
    for nx in nx_list:
        h  = L / nx
        dt = cfl * h / abs(v)

        sim = PureAdvection1D(
            c_velocity   = v,
            time_step    = dt,
            total_time   = T,
            spatial_step = nx
        )
        y_num = sim.run_simulation()
        x     = sim.x
        t_end = sim.Nt * sim.dt

        y_ex = gauss_exact(x, v, t_end)
        l2,_ = err(y_num, y_ex)

        h_list.append(h)
        e_list.append(l2)

        print(f"  Nx={nx:4d}, dx={h:.4f}, L2={l2:.6e}")

    h_arr = np.array(h_list)
    e_arr = np.array(e_list)

    # ----------- improved plot -----------
    plt.figure(figsize=(8,6), dpi=200)

    plt.loglog(h_arr, e_arr, "o-", linewidth=2, markersize=7,
               label="rms error (Gaussian pulse)")

    # slope=1 reference line
    ref_x = np.array([h_arr.min(), h_arr.max()])
    ref_y = e_arr[0] * (ref_x / h_arr[0])
    plt.loglog(ref_x, ref_y, "--", color="orange", linewidth=2,
               label="slope = 1")

    # axis labels (清晰、完整英文)
    plt.xlabel("step size (dx)")
    plt.ylabel("rms error (log scale)")

    plt.title("Convergence of first-order upwind (Gaussian pulse)")

    # major/minor ticks 更清楚
    plt.tick_params(axis="both", which="major", labelsize=12)
    plt.tick_params(axis="both", which="minor", labelsize=10)

    plt.grid(True, which="both", linestyle="--", alpha=0.4)

    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.show()


          


if __name__ == "__main__":
    all_plot()
    # resol_plot()  
    conv_plot()
