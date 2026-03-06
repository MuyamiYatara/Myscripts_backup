import numpy as np
import matplotlib.pyplot as plt
from itertools import islice
import os, re

# ===== 你需要自己维护这个映射表 =====
FILE_KIND_RULES = [
    (re.compile(r"seebeck_total_mu_.*eV\.dat$"), "seebeck"),
    (re.compile(r"rhotau_total_mu_.*eV\.dat$"), "rhotau"),
]

KIND_TITLE = {"seebeck": "Seebeck", "rhotau": r"\rho"}

# kind -> { col_index: (latex_label, unit_str, scale_factor) }
COLUMN_MAPS = {
    "seebeck": {
        # 例如：你说的 10 -> Szz，并且你以前乘 1e6（显示成 µV/K）
        10: (r"S_{zz}", r"\mu\mathrm{V}/\mathrm{K}", 1e6),
        6: (r"S_{yy}", r"\mu\mathrm{V}/\mathrm{K}", 1e6),
        2: (r"S_{xx}", r"\mu\mathrm{V}/\mathrm{K}", 1e6),
        # 你可以继续补：
        # 8: (r"S_{xx}", r"(\mu\mathrm{V}/\mathrm{K})", 1e6),
        # 9: (r"S_{yy}", r"(\mu\mathrm{V}/\mathrm{K})", 1e6),
    },
    "rhotau": {
        # 这里的单位/缩放你按你文件定义填
        # 举例：如果原始是 Ω·m，就可显示成 μΩ·cm：1 Ω·m = 1e8 μΩ·cm
        10: (r"\rho_{zz}", r"\mu\Omega\cdot\mathrm{cm}", 1e8),
        6: (r"\rho_{yy}", r"\mu\Omega\cdot\mathrm{cm}", 1e8),
        2: (r"\rho_{xx}", r"\mu\Omega\cdot\mathrm{cm}", 1e8),
        # 或者如果你的文件实际上是 rho*tau / tau 之类，也在这里写清楚
    }
}

# ===== 你需要自己维护这个映射表 =====




def get_col_meta(kind: str, ycol: int):
    m = COLUMN_MAPS.get(kind, {})
    if ycol not in m:
        raise KeyError(f"[{kind}] ycol={ycol} 不在 COLUMN_MAPS['{kind}'] 里，请补全映射。")
    return m[ycol]

def infer_kind_from_filename(file_name: str) -> str:
    base = os.path.basename(file_name)
    for pat, kind in FILE_KIND_RULES:
        if pat.search(base):
            return kind
    raise ValueError(f"无法从文件名识别 kind：{base}。请在 FILE_KIND_RULES 添加规则。")


def read_data(file_name: str):
    """
    读取类似：
      # T = 30 K
      <data block>
      # T = 40 K
      <data block>
    返回：dict[T] = np.array(block_data)
    """
    data = {}
    T = None
    block_data = []
    flag = 1

    with open(file_name, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if line.startswith("# T ="):
                if flag == 1:
                    flag = 0
                    T = float(line.split("=")[1].strip().split()[0])
                else:
                    data[T] = np.array(block_data)
                    T = float(line.split("=")[1].strip().split()[0])
                block_data = []
            elif line and not line.startswith("#"):
                block_data.append(list(map(float, line.split())))

        if T is not None:
            data[T] = np.array(block_data)

    return data



def plot_Y_vs_T_for_many_Btau(data_by_T: dict,
                             kind: str, 
                             ycol: int,
                             n_curves: int = 15,
                             stride: int = 1,
                             out_png: str = "Y-T.png",
                             title: str = "",
                             xlim=None,
                             ylim=None):
    """
    图1：横轴 T，曲线是不同 Btau（实际上是取固定 i 行：Brho[:, i, ycol]）
    - n_curves: 画多少条 Btau 曲线
    - stride: 每隔 stride 取一条（例如 stride=2 则 0,2,4,...）
    """
    # 按 T 排序
    Tlist = sorted(data_by_T.keys())
    blocks = [data_by_T[T] for T in Tlist]
    Brho = np.stack(blocks, axis=0)   # shape: (nT, nBtau, ncol)

    y_label_latex, y_unit, y_scale = get_col_meta(kind, ycol)

    plt.figure(figsize=(14, 8))
    nc = n_curves
    colors = plt.cm.viridis(np.linspace(0, 1, nc))

    # i 是 “Btau 行索引”
    # 你原来的写法是 i in count(start=0, step=st)，这里改成更直接的 range
    used = 0
    for k in range(n_curves):
        i = k * stride
        if i >= Brho.shape[1]:
            break

        y = Brho[:, i, ycol] * y_scale
        btau_val = Brho[0, i, 0]  # 默认第0列是 Btau（和你原脚本一致用法）
        plt.plot(Tlist, y, lw=3, color=colors[used], label=fr'$B\tau={btau_val:.3f}$')
        used += 1

    plt.xlabel('T (K)', fontsize=20)
    plt.ylabel(fr'${y_label_latex}\;({y_unit})$', fontsize=20)
    if title:
        plt.title(title, fontsize=24, pad=30)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(True)
    if xlim is not None:
        plt.xlim(*xlim)
    if ylim is not None:
        plt.ylim(*ylim)
    plt.legend(fontsize=14, loc='upper right', bbox_to_anchor=(1.3, 1))
    plt.subplots_adjust(right=0.78)
    plt.savefig(out_png, dpi=200)
    plt.close()


def plot_Y_vs_Btau_for_many_T(data_by_T: dict,
                              kind: str,
                              ycol: int,
                              T_take=(0, 28, 3),   # 等价你原来的 islice(...,0,28,3)
                              mode: str = "raw",   # "raw" 或 "norm0"
                              out_png: str = "Y-Btau.png",
                              title: str = "",
                              xlim=None,
                              ylim=None):
    """
    图2：横轴 Btau，曲线是不同 T
    mode:
      - "raw"  : 画 Y(Btau)
      - "norm0": 画 Y(Btau)/Y(0)
    """
    # 按 T 排序并抽样
    T_sorted = sorted(data_by_T.keys())
    start, stop, step = T_take
    T_pick = T_sorted[start:stop:step]

    y_label_latex, y_unit, y_scale = get_col_meta(kind, ycol)

    plt.figure(figsize=(14, 8))
    colors = plt.cm.viridis(np.linspace(0, 1, len(T_pick)))

    for i, T in enumerate(T_pick):
        arr = data_by_T[T]          # shape: (nBtau, ncol)
        x = arr[:, 0]               # 默认第0列是 Btau
        y = arr[:, ycol] * y_scale

        if mode == "raw":
            y_plot = y
            ylab = fr'${y_label_latex}\;({y_unit})$'
        elif mode == "norm0":
            y0 = y[0]
            # 避免除零
            y_plot = y / y0 if abs(y0) > 0 else np.full_like(y, np.nan)
            ylab = fr'${y_label_latex}(B\tau)/{y_label_latex}(0)$'
        else:
            raise ValueError("mode 只能是 'raw' 或 'norm0'")

        plt.plot(x, y_plot, lw=3, color=colors[i], label=f'T={T:.0f} K')

    plt.xlabel(r'$B(T)$', fontsize=20)
    plt.ylabel(ylab, fontsize=20)
    if title:
        plt.title(title, fontsize=24, pad=30)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(True)
    if xlim is not None:
        plt.xlim(*xlim)
    if ylim is not None:
        plt.ylim(*ylim)
    plt.legend(fontsize=14, loc='upper right', bbox_to_anchor=(1.3, 1))
    plt.subplots_adjust(right=0.78)
    plt.savefig(out_png, dpi=200)
    plt.close()


if __name__ == "__main__":
    # ====== 你只需要改这一小段配置就能换图 ======

    mu = 0.000
    file_name = f'./seebeck_total_mu_{mu:.3f}eV.dat'   # 或 rhotau_total_mu_...

    kind = infer_kind_from_filename(file_name)
    data = read_data(file_name)

    ycol = 10   # 你要画哪一列（例如 Szz / rhozz 对应的列）
    
    # 图1：Y(T)@不同Btau
    plot_Y_vs_T_for_many_Btau(
        data_by_T=data,
        kind=kind,
        ycol=ycol,
        out_png=f"{kind}-T.png",
        xlim=(0,300),
        title=fr'{kind}: col={ycol}, $\mu={mu*1000:.1f}\,\mathrm{{meV}}$'
    )

    # 图2：Y(Btau)@不同T（raw）
    plot_Y_vs_Btau_for_many_T(
        data_by_T=data,
        kind=kind,
        ycol=ycol,
        mode="raw",
        xlim=(0,1),
        out_png=f"{kind}-Btau.png",
    )

    # 图2：Y(Btau)/Y(0)（norm0）
    plot_Y_vs_Btau_for_many_T(
        data_by_T=data,
        kind=kind,
        ycol=ycol,
        mode="norm0",
        xlim=(0,1),
        ylim=(0.2,1.5),
        out_png=f"{kind}-Btau-norm0.png",
    )