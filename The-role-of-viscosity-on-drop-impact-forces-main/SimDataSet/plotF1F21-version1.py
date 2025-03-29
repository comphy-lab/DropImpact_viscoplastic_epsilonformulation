import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
import numpy as np

# ====================================================
# 配置 Matplotlib（使用LaTeX渲染）
# ====================================================
matplotlib.rcParams.update({
    "text.usetex": True,
    "font.size": 16,
    "axes.labelsize": 16,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 14,
    "text.latex.preamble": r"\usepackage{amsmath}"
})

# ====================================================
# 颜色方案（扩展支持更多Oh值）
# ====================================================
colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown', 'cyan', 'magenta']
color_map = {
    0.0025: colors[0],
    0.01: colors[1],
    0.025: colors[2],
    0.06: colors[3],
    0.1: colors[4],
    0.2: colors[5]
}

# ====================================================
# 数据读取（扩展Oh值范围）
# ====================================================
oh_values = [0.0025, 0.01, 0.06, 0.1, 0.2]  # 修正标点符号
dataframes = {}
for oh in oh_values:
    filename = f"Num_Oh{oh:.4f}.csv"  # 保持4位小数格式
    try:
        df = pd.read_csv(filename)
        if {"We", "F1", "F2"}.issubset(df.columns):
            dataframes[oh] = df
    except FileNotFoundError:
        print(f"Warning: {filename} not found. Skipping...")
        continue

# ====================================================
# 创建画布
# ====================================================
plt.figure(figsize=(8, 6))  # 增大画布尺寸以适应更大范围

# ====================================================
# 绘制数据（优化标记可见性）
# ====================================================
markers = ['o', 's', '^', 'D', 'v', 'p']  # 不同标记区分Oh值
lines_f1, labels_f1 = [], []
lines_f2, labels_f2 = [], []

for idx, (oh, df) in enumerate(dataframes.items()):
    # F1: 虚线 + 空心标记
    line_f1, = plt.plot(
        df["We"], df["F1"],
        marker=markers[idx], markersize=8,
        markeredgewidth=1, markerfacecolor='none',
        linestyle='--', linewidth=1.5,
        color=color_map[oh]
    )
    
    # F2: 实线 + 实心标记
    line_f2, = plt.plot(
        df["We"], df["F2"],
        marker=markers[idx], markersize=8,
        linestyle='-', linewidth=1.5,
        color=color_map[oh]
    )
    
    lines_f1.append(line_f1)
    labels_f1.append(rf'$Oh={oh}$')
    lines_f2.append(line_f2)
    labels_f2.append(rf'$Oh={oh}$')

# ====================================================
# 坐标轴设置（对数坐标优化）
# ====================================================
plt.xscale('log')  # 使用对数坐标
plt.yscale('linear')
plt.xlabel(r'$We$', fontsize=16)
plt.ylabel(r'$F/\rho V_0^2 D_0^2$', fontsize=16)

# 自定义刻度标签
plt.xlim(1, 100)
plt.xticks([1, 2, 5, 10, 20, 50, 100],
           [r'$1$', r'$2$', r'$5$', r'$10$', r'$20$', r'$50$', r'$100$'])

plt.ylim(0, 3)
plt.yticks(np.arange(0, 3.1, 0.5))

plt.grid(True, which='both', linestyle=':', alpha=0.5)

# ====================================================
# 图例优化（分栏显示）
# ====================================================
from matplotlib.legend import Legend

# 创建分栏图例
leg1 = Legend(plt.gca(), lines_f1, labels_f1,
              title=r'$F_1$',
              loc='upper left',
              bbox_to_anchor=(0.47, 0.98),
              frameon=False)

leg2 = Legend(plt.gca(), lines_f2, labels_f2,
              title=r'$F_2$',
              loc='upper left',
              bbox_to_anchor=(0.73, 0.98),
              frameon=False)

plt.gca().add_artist(leg1)
plt.gca().add_artist(leg2)

# ====================================================
# 保存输出
# ====================================================
plt.tight_layout()
plt.savefig("Extended_F1_F2_Plot.png", dpi=300, bbox_inches='tight')
plt.savefig("Extended_F1_F2_Plot.pdf", bbox_inches='tight')
# plt.show()