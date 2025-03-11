# 可视化结果
import matplotlib.pyplot as plt
import readfile

fsize = 30

#!--------------- plot a figure to compare the difference between the current parameter curve and the target curve ---------------!#
def compare_target_fit(X, Y_target, Y_fit, figname="res.png"):
    """
    X: X轴数据
    Y_target: 目标曲线的Y轴数据
    Y_fit: 当前参数下模型的Y轴数据
    figname: 保存图片的名称
    """

    plt.figure(figsize=(15, 9))

    plt.plot(X, Y_target, 'o', color='blue', linewidth=4 ,label='target' ,  alpha=0.5)
    plt.plot(X, Y_fit, 'r-', linewidth=4.5, label='fit')

    plt.xticks(fontsize=fsize)
    plt.yticks(fontsize=fsize)
    plt.title("Fit Result",fontsize=fsize)
    plt.xlabel("B(T)",fontsize=fsize)
    plt.ylabel("Y/(a.u.)",fontsize=fsize)
    
    plt.legend(fontsize=18, loc='upper right', bbox_to_anchor=(1.2, 1))  # 显示图例
    plt.subplots_adjust(right=0.85)
    plt.savefig(figname)
    plt.close()


if __name__ == "__main__":
    mat_name = "6V_101448_ScanB_I7-1Vxx5-3&8-10Vxy10-3&8-5ground2-11_hengya100-10k_RI1k__T1.9998K_Vbg6V_Vbias0V_ac1V.mat"
    B_and_target = readfile.read_matlab(mat_name, 'B', 'Rxx')
    print(len(B_and_target['B']), B_and_target['Rxx'])
    compare_target_fit(B_and_target['B'], B_and_target['Rxx'], B_and_target['Rxx'])

