import matplotlib.pyplot as plt

# 数据
y_values = [100, 10, 1, 0.1, 0.01, 0.001]
x_values = range(len(y_values))  # x轴为序号

# 绘图
plt.figure(figsize=(8, 5))
plt.plot(x_values, y_values, 'bo-', linewidth=2, markersize=8, label='Data')

# 设置对数刻度
plt.yscale('log')
plt.grid(True, which='both', linestyle='--', alpha=0.6)

# 添加公式注释
plt.text(0.5, 50, r'$Re = \frac{V_0 D_0}{\nu}$', fontsize=12, color='red')
plt.text(0.5, 5, r'$We = \frac{\rho V_0^2 D_0}{\gamma}$', fontsize=12, color='green')

# 标签和标题
plt.xlabel('Case Index', fontsize=12)
plt.ylabel('Value (log scale)', fontsize=12)
plt.title('Logarithmic Plot of Data with Re and We Formulas', fontsize=14)
plt.legend()

plt.show()