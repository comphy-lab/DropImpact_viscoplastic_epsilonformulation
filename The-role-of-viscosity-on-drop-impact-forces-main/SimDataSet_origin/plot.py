import matplotlib.pyplot as plt
import pandas as pd

# 指定 Oh 数值
oh_values = [0.0025, 0.01, 0.025, 0.06, 0.1, 0.2]

# 读取所有符合条件的 CSV 文件
dataframes = {}
for oh in oh_values:
    filename = f"Num_Oh{oh:.4f}.csv"
    try:
        df = pd.read_csv(filename)
        # 仅在包含 "We" 和 "F1" 时存入字典
        if {"We", "F2"}.issubset(df.columns):
            dataframes[oh] = df
        else:
            print(f"Warning: {filename} missing required columns, skipping.")
    except FileNotFoundError:
        print(f"Warning: {filename} not found.")

# 绘制 F2-We 图（x 轴对数）
plt.figure(figsize=(8, 6))
for oh, df in dataframes.items():
    if "We" in df.columns and "F2" in df.columns:
        plt.plot(df["We"], df["F2"], marker='s', linestyle='-', label=f"Oh={oh}")
plt.xscale("log")  # 设置 x 轴为对数坐标
plt.xlabel("We (log scale)")
plt.ylabel("F2")
plt.title("F2 vs We for Different Oh")
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.savefig("F2_We_plot.png", dpi=300)  # 保存图片
plt.show()
