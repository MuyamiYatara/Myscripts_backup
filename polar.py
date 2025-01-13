import matplotlib.pyplot as plt
import numpy as np

# 数据：角度（度数）和对应的半径（物理量）
data = np.loadtxt('combined_output.txt')
angles = data[:,0].astype(float)
radii = np.abs(data[:,10]).astype(float)

if np.any(np.isnan(radii)) or np.any(np.isinf(radii)):
    raise ValueError("radii contains NaN or inf values, which cannot be plotted.")

# 获取 radii 的最小值和最大值
r_min, r_max = np.min(radii), np.max(radii)

# 使用科学计数法设置yticks
yticks = np.linspace(0, r_max*0.1, 5).astype(float)  # 生成5个等间距的刻度
yticks_sci = ["{:.2e}".format(y) for y in yticks]  # 将yticks转换为科学计数法两位有效数字

# 提取科学计数法中的 10 的次方部分
exponent = int(np.floor(np.log10(r_max)))
scale_factor = 10 ** exponent

# 将 radii 数据归一化到10的次方
yticks_scaled = np.array(yticks / scale_factor, dtype=np.float64)
radii_scaled = np.array(radii / scale_factor, dtype=np.float64)


# 转换角度为弧度，因为matplotlib中的极坐标系统使用弧度
angles_rad = np.deg2rad(angles)
#print(angles_rad)

# 创建极坐标子图
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(10, 8))

# 绘制连线图
ax.plot(angles_rad, radii_scaled, linewidth=2)

# 设置极坐标网格（每30度一条线）
ax.set_theta_direction(-1)  # 使角度顺时针递增
ax.set_theta_offset(np.pi / 2)  # 将起点设为正上方（0度为顶部）

# 设置角度刻度标签
ax.set_xticks(np.deg2rad([0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]))
ax.set_xticklabels(['0°', '30°', '60°', '90°', '120°', '150°', '180°', '210°', '240°', '270°', '300°', '330°'])

# 隐藏极坐标图的默认径向刻度
ax.set_yticklabels([])

# 显示网格线
ax.grid(True)
# 设置径向范围
ax.set_ylim(0, r_max*1.1 / scale_factor)


# 创建一个次轴（普通轴）用于显示刻度
ax_y = fig.add_axes([0.1, 0.3, 0.05, 0.4])  # 在图形左边增加一个竖直轴
ax_y.spines['top'].set_color('none')
ax_y.spines['right'].set_color('none')
ax_y.spines['bottom'].set_color('none')
ax_y.spines['left'].set_position(('outward', 10))

# 设置次轴的刻度和标签
print(yticks_scaled)
ax_y.set_yticks(yticks_scaled)
ax_y.set_yticklabels(["{:.2f}".format(y) for y in yticks_scaled])
ax_y.set_ylabel(f"$\\times 10^{{{exponent}}}$", fontsize=14)
# 隐藏x轴的刻度线
ax_y.set_xticks([])


print("angles type:", angles.dtype)
print("radii type:", radii.dtype)
print("yticks type:", yticks.dtype)
print("yticks_scaled type:", yticks_scaled.dtype)
print("radii_scaled type:", radii_scaled.dtype)
print("angles_rad type:", angles_rad.dtype)
# 显示图形
plt.savefig("mr_polar.png")
