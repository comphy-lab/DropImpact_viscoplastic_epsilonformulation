import matplotlib.pyplot as plt
import pandas as pd

# 设置全局字体
plt.rcParams.update({'font.size': 18})

# 定义颜色方案 (与相图一致)
colors = {
    0.0025: 'royalblue',
    0.01: 'crimson',
    0.025: 'darkorange',
    0.06: 'green',
    0.1: 'purple',
    0.2: 'brown'
}

# 读取数据
oh_values = [0.0025, 0.01]
dataframes = {}
for oh in oh_values:
    filename = f"Num_Oh{oh:.4f}.csv"
    try:
        df = pd.read_csv(filename)
        if {"We", "F1", "F2"}.issubset(df.columns):
            dataframes[oh] = df
    except FileNotFoundError:
        print(f"Warning: {filename} not found.")

# 创建画布
plt.figure(figsize=(10, 8))

# 预存图例元素
f1_handles, f1_labels = [], []
f2_handles, f2_labels = [], []

# 绘制数据
for oh, df in dataframes.items():
    # 绘制F1
    line_f1, = plt.plot(df["We"], df["F1"], 
                        marker='s', markersize=8,
                        linestyle='-', linewidth=2,
                        color=colors[oh])
    f1_handles.append(line_f1)
    f1_labels.append(f"Oh={oh} (F1)")
    
    # 绘制F2
    line_f2, = plt.plot(df["We"], df["F2"], 
                        marker='o', markersize=8,
                        linestyle='--', linewidth=2,
                        color=colors[oh])
    f2_handles.append(line_f2)
    f2_labels.append(f"Oh={oh} (F2)")

# 设置坐标轴
# plt.xscale("log")
plt.xlim(12,30)
plt.ylim(0.3,1.2)
plt.xlabel("We", fontsize=18)
plt.ylabel("F", fontsize=18)
# plt.title("Viscoelastic Droplet Deformation", fontsize=20)
plt.grid(True, which="both", linestyle="--", alpha=0.6)

# 创建合并图例 (右侧内部)
leg = plt.legend(f1_handles + f2_handles, f1_labels + f2_labels,
                 loc='upper right',
                 title='Oh Values',
                 framealpha=0.9)

# 调整布局
plt.tight_layout()

plt.savefig("F1_F2_We_plot_inset.png", dpi=300, bbox_inches='tight')
plt.show()
