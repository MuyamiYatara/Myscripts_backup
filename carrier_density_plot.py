import numpy as np
from itertools import islice


def read_data_blocks(file_path):
    data_blocks = []  # 用于存储所有数据块
    current_block = []  # 用于存储当前数据块

    with open(file_path, 'r') as file:
        # 跳过前三行
        for line in islice(file, 3, None):
            stripped_line = line.strip()
            if not stripped_line:  # 空行表示一个数据块结束
                if current_block:
                    data_blocks.append(np.array(current_block, dtype=float))
                    current_block = []  # 清空当前块
            elif not stripped_line.startswith('#'):  # 跳过注释行
                current_block.append([float(x) for x in stripped_line.split()])
        
        # 添加最后一个块（如果存在）
        if current_block:
            data_blocks.append(np.array(current_block, dtype=float))
    
    return data_blocks

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

def plot_data_blocks(data_blocks, start, end):
    """
    绘制指定范围的数据块，每块数据用不同颜色标注。
    
    参数：
        data_blocks (list of np.ndarray): 数据块列表
        start (int): 要绘制的起始块编号
        end (int): 要绘制的结束块编号（包含）
    """
    # 创建颜色映射
    #colors = plt.cm.jet(np.linspace(0, 1, end - start + 1))
    num_colors = end - start + 1
    colors = []
    for i in range(num_colors):
        ratio = i / (num_colors - 1)  # 计算当前比例
        r = int(255 * (1 - ratio))   # 红色分量从 255 减少到 0
        g = int(255 * ratio)         # 绿色分量从 0 增加到 255
        b = 0                        # 蓝色分量保持为 0
        colors.append(f"#{r:02X}{g:02X}{b:02X}")  # 转换为十六进制颜色




    plt.figure(figsize= (18, 8))  # 设置图像大小
    # 获取当前坐标轴
    ax = plt.gca()

    # 使用 ScalarFormatter 设置科学计数法
    formatter = ScalarFormatter()
    formatter.set_scientific(True)
    formatter.set_powerlimits((0, 0))  # 始终启用科学计数法
    ax.yaxis.set_major_formatter(formatter)
    ax.yaxis.get_offset_text().set_fontsize(20)
    # 调整显示范围为 10^18
    for label in ax.get_yticklabels():
        label.set_text(f"{float(label.get_text()) / 1e18:.1f}e18")
    for i, block in enumerate(data_blocks[start:end + 1]):
        x = block[:, 0]  # 第一列作为 x
        y = block[:, 1]  # 第二列作为 y
        plt.plot(x, y, linewidth=3,label=f"mu = {(Minmu +  i*stepmu)*1000:.0f} meV", color=colors[i])
    
    plt.xlabel(r"T ($K$)",fontsize=21)
    plt.ylabel(r"n ($cm^{-3}$)",fontsize=21)
    plt.title(f"Carrier Density of Holes",fontsize=25) 
    #plt.title(f"Carrier Density of Electrons",fontsize=25)    
    #plt.title(f"carrier concentration of iband {((end+1)/Nummu):.0f}",fontsize=25)
    plt.xticks(fontsize=20)  # 设置x轴刻度字体大小
    plt.yticks(fontsize=20)
    plt.ylim(0.1*10**18,9*10**18)
    plt.legend(fontsize=18, loc='upper right', bbox_to_anchor=(1.3, 1))
    plt.subplots_adjust(right=0.80)
    plt.grid(True)
    plt.savefig(f'iband{((end+1)/Nummu):.0f}.png')





# 使用示例
NumT = 31
Minmu = -0.010
Maxmu = 0.010
Nummu = 11
iband = 5
stepmu = (Maxmu - Minmu)/(Nummu - 1.0)
file_path = 'carrier_concentration.dat'
data_blocks = read_data_blocks(file_path)

# 输出结果
#for i, block in enumerate(data_blocks):
#    print(f"数据块 {i}:")
#    print(block)

# 示例使用
# 假设 data_blocks 是之前读取的数组
# 从第 45 块到第 55 块绘图
plot_data_blocks(data_blocks, start=(iband-1)*Nummu, end=((iband)*Nummu-1))

