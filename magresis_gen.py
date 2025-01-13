from sklearn.metrics import silhouette_score, adjusted_rand_score, normalized_mutual_info_score
from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs
import numpy as np
import itertools
from joblib import Parallel, delayed
from tqdm import tqdm
import random
import matplotlib.pyplot as plt
import time
from scipy.optimize import fsolve


###------------------- 曲线生成与参数化 -------------------###
# 磁滞生成，alpha是磁滞饱和的速度，l越大磁滞面积越大
def M_BT_curve(B_array, T, alpha, l) :
    kb = 1.380649*10**-23
    mub = 9.28*10**-24
    Tc = 175.0 
    def equation(x):
        return x - np.tanh( mub*1.0/(kb*T)* B_v + (Tc / T) * x)
    
    lthB = len(B_array)
    MBT_array = []
    for i in B_array :
        B_v = i
        if i>0 :
            sols = fsolve(equation, 0.7)
        else :
            sols = -fsolve(equation, -0.7)
        MBT_array.append(sols[0]*np.tanh(alpha*B_v))
        


    dk = (MBT_array[l+1] - MBT_array[l])/(B_array[l+1] - B_array[l])
    b1 = MBT_array[l] - B_array[l]*dk
    b2 = MBT_array[len(MBT_array)-l-1] - B_array[lthB-l-1]*dk
    left_y =  dk*B_array[0] + b1
    right_y = dk*B_array[lthB-1] + b2

    MBT_loop = np.zeros(2*lthB)
    MBT_loop[0+lthB:lthB-l+lthB] = np.array(MBT_array[l:lthB])-(MBT_array[l]-left_y)
    MBT_loop[lthB-l+lthB:lthB+lthB] = np.array(B_array[lthB-l:lthB])*dk+b2
    MBT_loop[0:l] = np.array(B_array[0:l])*dk+b1
    MBT_loop[l:lthB] = np.array(MBT_array[0:lthB-l])-(MBT_array[lthB-1-l]-right_y)
    return MBT_loop


    


# 曲线定义
def mag_formula_two_band(sigma1, mu1, sigma2, mu2, sigmaA, B):


    #sigmaA_B = sigmaA*np.tanh(500*B)

    numerator = sigma1 + B**2*mu2**2*sigma1 + sigma2 + B**2*mu1**2*sigma2
    denominator = (1.0+B**2*mu2**2)*sigma1**2+ \
    2*sigma1*(sigma2 + B**2*mu1*mu2*sigma2 + B*mu1*(1.0+B**2*mu2**2)*sigmaA)+ \
    (1.0+B**2*mu1**2)*(sigma2**2 + 2*B*mu2*sigma2*sigmaA + (1.0+B**2*mu2**2)*sigmaA**2)

    return numerator/denominator

# 采样曲线的参数范围
curve_param_dictionary = { #每个参数都以 [最小值， 最大值， 个数] 的方式给定
    "B" : [-1, 1, 200], #单位：T
    "sigma1" : [1.0*10**4, 1.0*10**5, 10], #单位：S/cm 
    "mu1" : [30.0, 300.0, 10], #单位：(C*s)/kg, 为电子则mu为负值，为空穴则mu为正值
    "sigma2" : [1.0*10**4, 1.0*10**5, 10], #单位：S/cm
    "mu2" : [30.0, 300.0, 10], #单位：(C*s)/kg, 为电子则mu为负值，为空穴则mu为正值
    "sigmaA" : [100, 1000, 10], #单位：S/cm
}
paranamelist = list(curve_param_dictionary.keys())
# 定义每个参数对应的array
for pname, i in curve_param_dictionary.items() :
    globals()[pname + "_array"] = np.array([i[0] + j*(i[1]-i[0])/i[2] for j in range(i[2])])

# 生成总的参数空间
curve_param_array = np.array(list(itertools.product(sigma1_array, mu1_array, sigma2_array, mu2_array, sigmaA_array)))
total_curve_num = np.shape(curve_param_array)[0]


# 定义曲线生成函数，返回数组在磁场维度长度为两倍的len(B_array)
def generate_curve(start_idx, end_idx):
    #rhocurve_params = []
    #rhocurve_new_values = []
    rhocurve_values = []
    for i in range(start_idx, end_idx):
        lthB = len(B_array)
        y_1 = [mag_formula_two_band(curve_param_array[i, 0], curve_param_array[i, 1],\
                                   curve_param_array[i, 2], curve_param_array[i, 3], \
                                    curve_param_array[i, 4]*M_loop[index], j ) for index, j in enumerate(B_array)]
        
        y_2 = [mag_formula_two_band(curve_param_array[i, 0], curve_param_array[i, 1],\
                                   curve_param_array[i, 2], curve_param_array[i, 3], \
                                    curve_param_array[i, 4]*M_loop[index+lthB], j ) for index, j in enumerate(B_array)]
        
        y = np.hstack((np.array(y_1), np.array(y_2)))
        #归一化，使得他们的值都在[0,1]
        y_new = [ (j - min(y))/(max(y)- min(y)) for j in y]

        #print(y[101], B_array[101])
        #spline = UnivariateSpline(B_array[50:150], y[50:150], s=0, k=5)
        #y_new = spline(B_array)

        # 提取样条参数
        #rhocurve_params.append(spline.get_coeffs())  # 样条系数

        rhocurve_values.append(y_new)
        #rhocurve_new_values.append(y_new)
    return np.array(rhocurve_values)#, np.array(rhocurve_new_values), np.array(rhocurve_params)

# 使用 joblib 进行并行计算
n_jobs = 20 # 使用cpu数
curves_per_job = total_curve_num // n_jobs  # 每个 job 分配的曲线数量
remainder = total_curve_num % n_jobs  # 余数，即无法整除的曲线数
# 动态计算每个任务的起始和结束索引
job_ranges = []
start_idx = 0
for i in range(n_jobs):
    # 对前 'remainder' 个任务分配一个额外的曲线
    end_idx = start_idx + curves_per_job + (1 if i < remainder else 0)
    job_ranges.append((start_idx, end_idx))
    start_idx = end_idx

#生成磁滞回线
M_loop = M_BT_curve(B_array, 100.0, 5.0, 10)


# 使用 joblib 并行计算
results = Parallel(n_jobs=n_jobs)(
    delayed(generate_curve)(start, end)
    for start, end in tqdm(job_ranges)
)

# 汇总所有进程生成的曲线以及参数化后的结果
all_rhocurves = results#, all_rhocurves_new,  all_rhocurves_params = zip(*results) # 解包整体结果
all_rhocurves = np.concatenate(all_rhocurves, axis=0)
#all_rhocurves_params = np.concatenate(all_rhocurves_params, axis=0)
#all_rhocurves_new = np.concatenate(all_rhocurves_new, axis=0)


#random_integers = list(random.sample(range(100, 110), 3))
###############固定参数的调试###############
#i = 4435
#numnum = generate_curve(i, i+1)
#fsize = 30

# plt.figure(figsize=(15, 9)) 
# plt.plot(B_array, numnum[0,:],linewidth=3,linestyle='-', color='blue')
# plt.xlabel('B(T)',fontsize=fsize)  # x 轴标签
# plt.xticks(fontsize=fsize)
# plt.ylabel(r'$\rho_{xx}$($\Omega \cdot m$)',fontsize=fsize)  # y 轴标签
# plt.yticks(fontsize=fsize)
# plt.title(f'sigma1_{curve_param_array[i,0]:.1e}mu1_{curve_param_array[i,1]:.1e}'+\
#           f'sigma2_{curve_param_array[i,2]:.1e}mu2_{curve_param_array[i,3]:.1e}sigmaA_{curve_param_array[i,4]:.1e}', fontsize=24, pad=30)  # 图像标题
# #plt.legend(fontsize=25, loc='upper right', bbox_to_anchor=(1.35, 1))  # 显示图例
# plt.subplots_adjust(right=0.9)
# plt.tick_params(pad=8)
# #plt.xlim(0,0.5)
# #plt.ylim(0,40)
# # 显示图像
# plt.grid(True)  # 添加网格
# plt.savefig(f"{i}.png")
###############固定参数的调试###############



# for i in random_integers :
#     fsize = 30
#     print(i)

#     print("params", all_rhocurves_params[i,:])
    
#     plt.figure(figsize=(15, 9)) 
#     plt.plot(B_array, all_rhocurves[i,:],linewidth=3,linestyle='-', color='blue',label="curve")
#     plt.plot(B_array, all_rhocurves_new[i,:],linewidth=3,linestyle='-', color='red',label="fit")
#     plt.xlabel('B(T)',fontsize=fsize)  # x 轴标签
#     plt.xticks(fontsize=fsize)
#     plt.ylabel(r'$\rho_{xx}$($\Omega \cdot m$)',fontsize=fsize)  # y 轴标签
#     plt.yticks(fontsize=fsize)
#     plt.title(f'sigma1_{curve_param_array[i,0]:.1e}mu1_{curve_param_array[i,1]:.1e}'+\
#               f'sigma2_{curve_param_array[i,2]:.1e}mu2_{curve_param_array[i,3]:.1e}sigmaA_{curve_param_array[i,4]:.1e}', fontsize=24, pad=30)  # 图像标题
#     plt.legend(fontsize=18, loc='upper right', bbox_to_anchor=(1.2, 1))  # 显示图例
#     plt.subplots_adjust(right=0.85)
#     plt.tick_params(pad=8)
#     #plt.xlim(0,0.5)
#     #plt.ylim(0,40)
#     # 显示图像
#     plt.grid(True)  # 添加网格
#     plt.savefig(f"{i}_gen.png")



# 保存最终结果
#np.save('generated_curves.npy', all_rhocurves)



###------------------- 曲线生成与参数化 -------------------###


###------------------- 聚类 -------------------###

# 聚类
def fit_kmeans(par, n_clusters) :
    kmeans = KMeans(n_clusters=n_clusters ,random_state=0).fit(par)
    labels = kmeans.labels_
    silhouette_avg = silhouette_score(par, labels)
    return n_clusters, labels, silhouette_avg

n_clusters_list = [13, 12, 13]
kmeans_results = Parallel(n_jobs=3)(delayed(fit_kmeans)(all_rhocurves, n_clusters) for n_clusters in tqdm(n_clusters_list))
n_list, label_list, avg_list  = zip(*kmeans_results)
print(n_list, avg_list)


labels = label_list[0]
unique_labels = np.unique(labels)  # 获取所有唯一的聚类标签
all_rhocurves_clusters = []
all_rhocurves_clusters_index = []
for label in unique_labels:
    indices = np.where(labels == label)[0]  # 获取该类的所有索引
    cluster_data = all_rhocurves[indices]  # 提取该类对应的数据
    all_rhocurves_clusters.append(cluster_data)
    all_rhocurves_clusters_index.append(indices)

    


# # 创建一个字典来存储每个聚类的索引
# cluster_indices = defaultdict(list)

# for idx, label in enumerate(labels):
#     cluster_indices[label].append(idx)

# 内部评价
silhouette_avg = silhouette_score(all_rhocurves, labels)
print("轮廓系数:", silhouette_avg)
print("各类指标:", all_rhocurves_clusters_index)

# 第k类
for k in range(13) :
    for i,index in enumerate(all_rhocurves_clusters_index[k][0:11]) :
        fsize = 30
        
        plt.figure(figsize=(15, 9)) 
        plt.plot(B_array, all_rhocurves_clusters[k][i,0:len(B_array)],linewidth=3,linestyle='-', color='blue',label="curve")
        plt.plot(B_array, all_rhocurves_clusters[k][i,len(B_array):len(B_array)*2],linewidth=3,linestyle='-', color='blue',label="curve")
        plt.xlabel('B(T)',fontsize=fsize)  # x 轴标签
        plt.xticks(fontsize=fsize)
        plt.ylabel(r'$\rho_{xx}$($\Omega \cdot m$)',fontsize=fsize)  # y 轴标签
        plt.yticks(fontsize=fsize)
        plt.title(f'sigma1_{curve_param_array[index,0]:.1e}mu1_{curve_param_array[index,1]:.1e}'+\
                f'sigma2_{curve_param_array[index,2]:.1e}mu2_{curve_param_array[index,3]:.1e}sigmaA_{curve_param_array[index,4]:.1e}', fontsize=24, pad=30)  # 图像标题
        plt.legend(fontsize=18, loc='upper right', bbox_to_anchor=(1.2, 1))  # 显示图例
        plt.subplots_adjust(right=0.85)
        plt.tick_params(pad=8)
        #plt.xlim(0,0.5)
        #plt.ylim(0,40)
        # 显示图像
        plt.grid(True)  # 添加网格
        plt.savefig(f"k={k}_{index}.png")

###------------------- 聚类 -------------------###






