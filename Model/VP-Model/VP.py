import numpy as np
import matplotlib.pyplot as plt

# 定义X范围（log空间）
x = np.logspace(-2, 1, 500)  # 从0.001到10

# 定义不同的J值
J_values = [0.1,1,10]

# 绘图
plt.figure(figsize=(8, 6))
for J in J_values:
    y = (J / (x + 0.001*J) + 0.0035) * x
    y_ref = J + 0.0035 * x  # 参考虚线
    plt.plot(x, y, label=f'Y (J={J})')
    plt.plot(x, y_ref, '--', label=f'Y = ({J} + 0.001) * X', alpha=0.7)

# 设置log坐标轴
plt.xscale('log')
plt.yscale('log')

# 添加图例和标签
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Y = (J / (X + 0.01 * J)) * X in log-log space')
plt.legend()
plt.grid(True, which="both", ls="--", lw=0.5)

plt.tight_layout()
plt.savefig("VP.png",dpi=300)
