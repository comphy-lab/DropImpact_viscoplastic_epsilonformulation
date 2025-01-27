import os
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
We = 20
Js=[0,0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
Ohs=map(lambda x: f"{x:.4f}", np.logspace(-3, 0, 16))
MAXlevels = [9]

# 结果文件所在的目录
model_name=f"Results_Running_We{We}"
currentdir=os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.join(currentdir, model_name) # 替换为你的结果文件所在目录

# 创建一个列表来存储所有结果
summary_data = []

# 遍历所有参数组合

for MAXlevel in MAXlevels:
    for Oh in Ohs:
        for J in Js:
            # 构建文件名
            filename = f"Bo0.0-We{We}-J{J}-Oh{Oh}-MAXlevel{MAXlevel}-epsilon0.01.csv"
            file_path = os.path.join(results_dir, filename)
            
            # 检查文件是否存在
            if not os.path.isfile(file_path):
                print(f"{file_path} 不存在，跳过。")
                continue
            
            # 读取数据文件
            try:
                data = pd.read_csv(file_path)
                data.fillna(0, inplace=True)
            except Exception as e:
                print(f"读取 {filename} 时出错：{e}")
                continue
            
            # 检查必要的列是否存在
            required_columns = ["t", "Zmin"]
            if not all(col in data.columns for col in required_columns):
                print(f"{filename} 缺少必要的列，跳过。")
                continue
            
            min_Zmin = data["Zmin"].min()
            Zmin_last = data["Zmin"].iloc[-1]-min_Zmin

            # 确定 phase 的值
            if Zmin_last < 0.02:
                phase = 1
            else:            
                phase = 0

            # 将结果添加到 summary_data 列表
            summary_data.append({
                "Oh": Oh,
                "J": J,
                "Phase": phase
            })

# 将结果转换为 DataFrame
summary_df = pd.DataFrame(summary_data)

# 按需要保存结果，例如保存为 CSV 文件
summary_df.to_csv(os.path.join(currentdir, f"We{We}_Phase.csv"), index=False)

import matplotlib.pyplot as plt

# 绘制 phase 相图 (J vs Oh)，使用 log 坐标
color_map = {1: 'blue', 0: 'red'}
plt.figure(figsize=(10, 10))
scatter = plt.scatter(
    summary_df['Oh'],
    summary_df['J'],
    c=summary_df['Phase'].map(color_map),  # 映射颜色
    s=100,
    edgecolor='k'
)

# plt.xscale("log")

# 添加标签和标题
plt.xlabel("Oh")
plt.ylabel("J")  
plt.title(f"We{We}-Phase Diagram in J-Oh Space")

# 显示图像
plt.grid(which="both", linestyle="--", alpha=0.5)
plt.savefig(os.path.join(currentdir, f"We{We}_Phase.png"), dpi=300, bbox_inches='tight')
plt.show()
