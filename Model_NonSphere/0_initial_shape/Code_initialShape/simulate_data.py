import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# 设定文件夹路径
folder_path = "simulation"

# 设定 J 和 R_y 参数
Js = [0.01, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0]
R_ys = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# 结果存储
results = []

# 遍历所有组合
for J in Js:
    for R_y in R_ys:
        filename = f"Bo0.0-We10-J{J}-R_y{R_y}-Oh0.01-MAXlevel9-epsilon0.001.csv"
        file_path = os.path.join(folder_path, filename)
        
        if not os.path.exists(file_path):
            continue
        
        # 读取数据
        df = pd.read_csv(file_path)
        
        # 确保数据包含所需列
        if {'ke', 'Rmax'}.issubset(df.columns):
            ke_max = df['ke'].max()
            Rmax_max = df['Rmax'].max()
            Rmax_min = df['Rmax'].min()
            Delta = Rmax_max - Rmax_min
            phase = 0 if Delta < 0.01 else 1
            results.append([J, R_y, ke_max, Delta, Rmax_max, Rmax_min, phase])

# 转换为 DataFrame
output_df = pd.DataFrame(results, columns=["J", "R_y", "ke_max", "Delta", "Rmax_max", "Rmax_min", "phase"])

# 导出 CSV
output_csv_path = "simulation_results.csv"
output_df.to_csv(output_csv_path, index=False)

print(f"结果已保存至 {output_csv_path}")

# 绘制 phase 图
plt.figure(figsize=(12, 10))  # 加宽画布容纳右侧图例

# 分离不同相位的数据
stable_df = output_df[output_df["phase"] == 0]
unstable_df = output_df[output_df["phase"] == 1]

# 分别绘制不同相位（调整标记大小）
plt.scatter(stable_df["R_y"], stable_df["J"], 
            c='royalblue', edgecolors='k', s=120,
            label='Stable', zorder=3)
plt.scatter(unstable_df["R_y"], unstable_df["J"], 
            c='crimson', edgecolors='k', s=120,
            label='Unstable', zorder=3)

plt.xscale("log")
# plt.yscale("log")
plt.xlabel(r"$R_{\theta}/R_0$", fontsize=18)
plt.ylabel(r"$J$", fontsize=18)
plt.title("Stability of the elliptical viscoplastic droplet", fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

# 绘制理论曲线（修正公式显示）
R_y_curve = np.logspace(-1, 1, 400)
J_curve = 0.6 * np.abs(R_y_curve**2 - 1) / R_y_curve
plt.plot(R_y_curve, J_curve, color='black', linewidth=2.5, 
         linestyle='--', label=r'Criteria: $J=0.6\left|\frac{R_\theta^2 - 1}{R_\theta}\right|$')

# 设置图例（右侧垂直排列）
plt.legend(
    loc='upper center',
    bbox_to_anchor=(0.5, -0.12),  # 将图例放在图表下方
    ncol=3,  # 三列布局
    frameon=True,
    fontsize=18,
    framealpha=0.9,
    borderpad=1,
    borderaxespad=1
)

# 设置网格和布局
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout(rect=[0, 0, 0.85, 1])  # 右侧留出15%空间

plt.savefig("phase_diagram.png", dpi=300, bbox_inches='tight')
plt.show()