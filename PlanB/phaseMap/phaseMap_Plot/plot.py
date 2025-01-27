import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob

# 定义参数
Js = np.array(["0", "0.005", "0.01", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"])
Ohs = np.logspace(-3, 0, 16)
We=20
# 创建结果存储的DataFrame
results = []

# 遍历文件夹中的所有文件
folder_path = "mapWe20ResultsLevel9"
for J in Js:
    for Oh in Ohs:
        # 构建文件名模式
        file_path = f"{folder_path}/Bo0.0-We{We}-J{J}-Oh{Oh:.4f}-MAXlevel9-epsilon0.01.csv"
        # 读取CSV文件
        df = pd.read_csv(file_path)
        
        # 计算maxvc
        maxvc = df['vc'].max()
        
        # 计算phase
        Zmin = df['Zmin']
        phase = 1 if (Zmin.iloc[-1] - Zmin.min()) > 0.02 else 0
        
        # 存储结果
        results.append({
            'J': J,
            'Oh': Oh,
            'maxvc': maxvc,
            'phase': phase
        })

# 创建结果DataFrame并保存为CSV
results_df = pd.DataFrame(results)
results_df.to_csv('J-Oh-maxvc-phase.csv', index=False)

# 创建图形
plt.figure(figsize=(12, 12))
# J-Oh-maxvc相图
scatter = plt.scatter(results_df['Oh'], results_df['J'], 
                     c=results_df['maxvc'], cmap='viridis')
plt.colorbar(scatter, label='maxvc')
plt.ylabel('J')
plt.xlabel('Oh')
plt.xscale("log")
plt.title('J-Oh-maxvc Phase Diagram')
plt.tight_layout()
plt.savefig("J-Oh-maxvc.png",dpi=300)
plt.close()

# J-Oh-phase散点图
plt.figure(figsize=(12, 12))
colors = ['blue' if p == 0 else 'red' for p in results_df['phase']]
plt.scatter(results_df['Oh'],results_df['J'], c=colors)
plt.ylabel('J')
plt.xlabel('Oh')
plt.xscale("log")
plt.title('J-Oh-phase Scatter Plot')
plt.tight_layout()
plt.savefig("J-Oh-phase.png",dpi=300)
plt.close()