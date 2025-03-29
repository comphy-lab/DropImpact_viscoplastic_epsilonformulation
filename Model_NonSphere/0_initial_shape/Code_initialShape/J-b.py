import numpy as np
import matplotlib.pyplot as plt
import os

# 生成数据
b = np.logspace(-1, 1, 400)  # 0.1到10对数间隔
J = np.abs(b**2 - 1) / b

# 创建画布
plt.figure(figsize=(10, 6))

# 绘制曲线
main_line = plt.plot(b, J, color='darkblue', linewidth=3, zorder=3)

# 填充区域
plt.fill_between(b, J, 0, color='gold', alpha=0.2, label='Stable Region')  # 下方淡黄
plt.fill_between(b, J, 4, color='lightblue', alpha=0.3, label='Unstable Region')  # 上方淡蓝

# 坐标轴设置
plt.xscale('log')
plt.xlim(0.2, 5)
plt.ylim(0, 4)
plt.xlabel(r'$R_{\theta}$', fontsize=24, labelpad=10)
plt.ylabel('J', fontsize=24, rotation=0, labelpad=10)

# 刻度设置
plt.xticks(fontsize=18)
plt.yticks(np.arange(0, 5, 1), fontsize=18)
plt.grid(True, which='major', linestyle='--', linewidth=1, alpha=0.7)

# 保存图片
script_dir = os.path.dirname(os.path.abspath(__file__))
save_path = os.path.join(script_dir, "stability_plot.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')  # 保存为PNG

# 显示图形
plt.show()
