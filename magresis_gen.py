from sklearn.metrics import silhouette_score, adjusted_rand_score, normalized_mutual_info_score
from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs
from openTSNE import TSNE
import numpy as np
import itertools
from joblib import Parallel, delayed
from tqdm import tqdm
import random
import matplotlib.pyplot as plt
import time
import sys
from scipy.optimize import fsolve

kb = 1.380649*10**-23 #单位：J/K
hbar = 1.0545718*10**-34 #单位：J*s
temperature = 100.0  #单位：K
electron = 1.6*10**-19 #单位：C
trans_pt = 5 #控制磁滞曲线的面积，从哪个B点开始转变磁滞
gnumber = 12 #类数目
alpha = 50000

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
    for d, i in enumerate(B_array) :
        
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

def mag_formula_two_band(n1, mu1, n2, mu2, sigmaA, B):

    sigma1 = electron*n1*np.abs(mu1)
    sigma2 = electron*n2*np.abs(mu2)

    sigma_matrix_1 = np.array([[sigma1/(1+mu1**2*B**2), sigma1*mu1*B/(1+mu1**2*B**2)],\
                              [-sigma1*mu1*B/(1+mu1**2*B**2), sigma1/(1+mu1**2*B**2)]])
    
    sigma_matrix_2 = np.array([[sigma2/(1+mu2**2*B**2), sigma2*mu2*B/(1+mu2**2*B**2)],\
                              [-sigma2*mu2*B/(1+mu2**2*B**2), sigma2/(1+mu2**2*B**2)]])
    
    sigma_matrix_A = np.array([[0.0, sigmaA],\
                              [-sigmaA, 0.0]])

    sigma_matrix_total= sigma_matrix_1 + sigma_matrix_2 + sigma_matrix_A

    #sigma_inv = (1.0/np.linalg.det(sigma_matrix_total))*np.transpose( sigma_matrix_total)#
    sigma_inv = np.linalg.inv(sigma_matrix_total)
    return sigma_inv[0,0]


def mag_formula_three_band(sigma1, mu1, sigma2, mu2, sigma3, mu3, sigmaA, B):
    sigma_matrix_1 = np.array([[sigma1/(1+mu1**2*B**2), sigma1*mu1*B/(1+mu1**2*B**2)],\
                              [-sigma1*mu1*B/(1+mu1**2*B**2), sigma1/(1+mu1**2*B**2)]])
    
    sigma_matrix_2 = np.array([[sigma2/(1+mu2**2*B**2), sigma2*mu2*B/(1+mu2**2*B**2)],\
                              [-sigma2*mu2*B/(1+mu2**2*B**2), sigma2/(1+mu2**2*B**2)]])

    sigma_matrix_3 = np.array([[sigma3/(1+mu3**2*B**2), sigma3*mu3*B/(1+mu3**2*B**2)],\
                              [-sigma3*mu3*B/(1+mu3**2*B**2), sigma3/(1+mu3**2*B**2)]])
    sigma_matrix_A = np.array([[0.0, sigmaA],\
                              [-sigmaA, 0.0]])

    sigma_matrix_total= sigma_matrix_1 + sigma_matrix_2 + sigma_matrix_3 + sigma_matrix_A
    sigma_inv = np.linalg.inv(sigma_matrix_total)
    return sigma_inv[0,0]


# def magthermo_formula_two_band(sigma1, mu1, sigma2, mu2, sigmaA, B) :
#     sigma_matrix_1 = np.array([[sigma1/(1+mu1**2*B**2), sigma1*mu1*B/(1+mu1**2*B**2)],\
#                               [-sigma1*mu1*B/(1+mu1**2*B**2), sigma1/(1+mu1**2*B**2)]])
    
#     sigma_matrix_2 = np.array([[sigma2/(1+mu2**2*B**2), sigma2*mu2*B/(1+mu2**2*B**2)],\
#                               [-sigma2*mu2*B/(1+mu2**2*B**2), sigma2/(1+mu2**2*B**2)]])

                        

# 采样曲线的参数范围
curve_param_dictionary = { #每个参数都以 [最小值， 最大值， 个数] 的方式给定
    #!# 请确保B的数值个数为奇数，否则将会带来奇异图像以及特征计算错误
    "B" : [-0.5, 0.5, 201], #单位：T
    #!# 请确保B的数值个数为奇数，否则将会带来奇异图像以及特征计算错误

    "n1" : [1.0*10**25, 5.0*10**25, 20], #单位：m^-3 
    "mu1" : [-11, 10, 20], #单位：(C*s)/kg, 为电子则mu为负值，为空穴则mu为正值
    "n2" : [1.0*10**25, 5.0*10**25, 20], #单位：m^-3
    "mu2" : [-11, 10, 20], #单位：(C*s)/kg, 为电子则mu为负值，为空穴则mu为正值
    #"sigma3" : [1.0*10**5, 5.0*10**5, 2], #单位：S/m
    #"mu3" : [-300.0, 300.0, 30], #单位：(C*s)/kg, 为电子则mu为负值，为空穴则mu为正值
    "sigmaA" : [1.0*10**7, 2000, 1], #单位：S/m
    #"tau": [1.0*10**-13, 10.0*10**-13, 2],
    #"seebeckA" :[1.0, 3.0, 2]
}
paranamelist = list(curve_param_dictionary.keys())
# 定义每个参数对应的array
for pname, i in curve_param_dictionary.items() :
    globals()[pname + "_array"] = np.array([i[0] + j*(i[1]-i[0])/i[2] for j in range(i[2])])

# 生成总的参数空间

curve_param_array = np.array(list(itertools.product(n1_array, mu1_array, \
                                                    n2_array, mu2_array, sigmaA_array ))) #@@@ two bands

#curve_param_array = np.array(list(itertools.product(sigma1_array, mu1_array, sigma2_array, \
#                                                    mu2_array, sigma3_array, mu3_array,sigmaA_array ))) #@@@ three bands

total_curve_num = np.shape(curve_param_array)[0]

# 定义曲线的特征函数，包括曲线数值的归一化，以及从磁滞转变点开始的斜率、二阶导数、磁滞转变点的数值
def characterize(yvals, B_lenth) :

    y1_front = yvals[0:int((B_lenth+1)/2)+trans_pt]
    y1_back = yvals[int((B_lenth+1)/2)+trans_pt:int(B_lenth)]
    y2_front = yvals[int(B_lenth):int((B_lenth+1)*3/2)-trans_pt] 
    y2_back = yvals[int((B_lenth+1)*3/2)-trans_pt:int(B_lenth*2)] #提取正负两条曲线的前后半段分别提取特征，这样每个部分都是连续的

    #一阶导数以及断点的值
    y1_front_de1 = [(y1_front[n+1] - y1_front[n]) for n in range(len(y1_front)-1)]
    y2_front_de1 = [(y2_front[n+1] - y2_front[n]) for n in range(len(y2_front)-1)]
    y1_front_de1.append(y1_front[len(y1_front)-1])
    y2_front_de1.append(y2_front[len(y2_front)-1]) 

    y1_back_de1 = [(y1_back[n] - y1_back[n-1]) for n in range(1,len(y1_back))]
    y2_back_de1 = [(y2_back[n] - y2_back[n-1]) for n in range(1,len(y2_back))]
    y1_back_de1.insert(0, y1_back[0])
    y2_back_de1.insert(0, y2_back[0]) 

    y1_front_de2 = [(y1_front[n+1] - 2*y1_front[n] +y1_front[n-1]) for n in range(1, len(y1_front)-1)]
    y2_front_de2 = [(y2_front[n+1] - 2*y2_front[n] +y2_front[n-1]) for n in range(1, len(y2_front)-1)]
    y1_back_de2 = [(y1_back[n+1] - 2*y1_back[n] +y1_back[n-1]) for n in range(1, len(y1_back)-1)]
    y2_back_de2 = [(y2_back[n+1] - 2*y2_back[n] +y2_back[n-1]) for n in range(1, len(y2_back)-1)]

    y1_front_chara =np.hstack((np.array(y1_front_de1), np.array(y1_front_de2)))
    y2_front_chara =np.hstack((np.array(y2_front_de1), np.array(y2_front_de2)))
    y1_back_chara =np.hstack((np.array(y1_back_de1), np.array(y1_back_de2)))
    y2_back_chara =np.hstack((np.array(y2_back_de1), np.array(y2_back_de2)))



    y_characters = np.hstack((np.array(y1_front_chara),np.array(y1_back_chara),np.array(y2_front_chara),np.array(y2_back_chara)))
    #print (np.shape(y_characters))
    return y_characters


# 定义曲线生成函数，返回数组在磁场维度长度为两倍的len(B_array)
def generate_curve(start_idx, end_idx):
    #rhocurve_params = []
    #rhocurve_new_values = []
    rhocurve_values = []
    rhocurve_characters = []
    for i in range(start_idx, end_idx):
        lthB = len(B_array)
        y_1 = [mag_formula_two_band(curve_param_array[i, 0], curve_param_array[i, 1],\
                                    curve_param_array[i, 2], curve_param_array[i, 3], \
                                     curve_param_array[i, 4]*M_loop[index], j ) for index, j in enumerate(B_array)] #@@@ two bands
        
        y_2 = [mag_formula_two_band(curve_param_array[i, 0], curve_param_array[i, 1],\
                                    curve_param_array[i, 2], curve_param_array[i, 3], \
                                     curve_param_array[i, 4]*M_loop[index+lthB], j ) for index, j in enumerate(B_array)] #@@@ two bands
        
        #y_1 = [mag_formula_three_band(curve_param_array[i, 0], curve_param_array[i, 1],\
        #                           curve_param_array[i, 2], curve_param_array[i, 3], \
        #                           curve_param_array[i, 4], curve_param_array[i, 5], \
        #                            curve_param_array[i, 6]*M_loop[index], j ) for index, j in enumerate(B_array)] #@@@ three bands
        
        #y_2 = [mag_formula_three_band(curve_param_array[i, 0], curve_param_array[i, 1],\
        #                           curve_param_array[i, 2], curve_param_array[i, 3], \
        #                           curve_param_array[i, 4], curve_param_array[i, 5], \
        #                            curve_param_array[i, 6]*M_loop[index+lthB], j ) for index, j in enumerate(B_array)] #@@@ three bands
        y = np.hstack((np.array(y_1), np.array(y_2)))
        #y_new = [((j - np.min(y))/(np.max(y)- np.min(y))) for j in y]
        #归一化，使得他们的值都在[0,1]
        y_new = []
        y_min = np.min(y)
        y_max = np.max(y)
        y_range = y_max - y_min

        # 使用y_range直接判断并执行
        #y_new = [(j - y_min) / y_range if y_range > 1e-7 else (j / y_max) for j in y]
        # 带反常的不会出现常数磁阻曲线，所以不判断了
        y_new = [(j - y_min) / y_range for j in y]
        dy = characterize(yvals=y_new, B_lenth=lthB)


        # 提取样条参数
        #rhocurve_params.append(spline.get_coeffs())  # 样条系数
        rhocurve_values.append(y_new)
        rhocurve_characters.append(y_new)
        #rhocurve_new_values.append(y_new)
    return np.array(rhocurve_values), np.array(rhocurve_characters)#, np.array(rhocurve_new_values), np.array(rhocurve_params)

# 使用 joblib 进行并行计算
n_jobs =50 # 使用cpu数
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
M_loop = M_BT_curve(B_array, temperature, alpha, trans_pt)


# 使用 joblib 并行计算
results = Parallel(n_jobs=n_jobs)(
    delayed(generate_curve)(start, end)
    for start, end in tqdm(job_ranges)
)

# 汇总所有进程生成的曲线以及参数化后的结果
all_rhocurves, all_rhocurves_characters = zip(*results)#, all_rhocurves_new,  all_rhocurves_params = zip(*results) # 解包整体结果
all_rhocurves = np.concatenate(all_rhocurves, axis=0)
all_rhocurves_characters = np.concatenate(all_rhocurves_characters, axis=0)
#-------------##-------------##-------------##-------------##-------------##-------------#

# fsize=30
# index=0
# g=0
# print(f"{all_rhocurves[index,0:len(B_array)]}", "###", all_rhocurves[0,205], "###", M_loop[int((len(B_array)+1)/2)+trans_pt-1: int((len(B_array)+1)/2)+trans_pt+1])
# #all_rhocurves = np.array([[1*i/i for i in range(1,401)]])
# # all_rhocurves = all_rhocurves[0]
# plt.figure(figsize=(15, 9)) 

# #all_rhocurves = np.array(all_rhocurves)
# plt.plot(B_array, all_rhocurves[index,0:len(B_array)],linewidth=3,linestyle='-', color='blue',label="+")
# plt.plot(B_array, all_rhocurves[index,len(B_array):len(B_array)*2],linewidth=3,linestyle='-', color='red',label="-")
# plt.xlabel('B(T)',fontsize=fsize)  # x 轴标签
# plt.xticks(fontsize=fsize)
# plt.ylabel(r'$\rho_{xx}$($a.u.$)',fontsize=fsize)  # y 轴标签
# plt.yticks(fontsize=fsize)
# plt.title(f'n1_{curve_param_array[index,0]:.3e}mu1_{curve_param_array[index,1]:.3e}'+\
#                 f'n2_{curve_param_array[index,2]:.3e}mu2_{curve_param_array[index,3]:.3e}'+\
#                     f'sigmaA_{curve_param_array[index,4]:.1e}', fontsize=20, pad=30)  # 图像标题
# plt.legend(fontsize=18, loc='upper right', bbox_to_anchor=(1.2, 1))  # 显示图例
# plt.subplots_adjust(right=0.85)
# plt.tick_params(pad=8)
# #plt.xlim(0,0.5)
# plt.ylim(0,1)
# # 显示图像
# plt.grid(True)  # 添加网格
# plt.savefig(f"compare.png")
# sys.exit()
#-------------##-------------##-------------##-------------##-------------##-------------#
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
    #silhouette_avg = silhouette_score(par, labels)
    return n_clusters, labels#, silhouette_avg

n_clusters_list = [gnumber]
print(np.array(all_rhocurves).shape)
kmeans_results = Parallel(n_jobs=1)(delayed(fit_kmeans)(all_rhocurves_characters, n_clusters) for n_clusters in tqdm(n_clusters_list))
n_list, label_list  = zip(*kmeans_results)
print(n_list)#, avg_list)


labels = label_list[0]
unique_labels = np.unique(labels)  # 获取所有唯一的聚类标签
all_rhocurves_clusters = []
all_rhocurves_clusters_index = []
for label in unique_labels:
    indices = np.where(labels == label)[0]  # 获取该类的所有索引
    cluster_data = all_rhocurves[indices]  # 提取该类对应的数据
    all_rhocurves_clusters.append(cluster_data)
    all_rhocurves_clusters_index.append(indices)

    
# 使用tsne对数据进行降维
tsne = TSNE(n_components=2, 
            random_state=20, 
            perplexity=5, 
            n_iter=40, 
            metric="euclidean", 
            n_jobs=1)  # 默认初始化方法
curves_tsne = tsne.fit(all_rhocurves_characters)

# 绘制tsne降维后的散点图


# 绘制tsne结果

plt.figure(figsize=(12, 5),dpi=150)
fsize = 30
scatter = plt.scatter(curves_tsne[:, 0], curves_tsne[:, 1], c=labels, cmap='tab20', s=10)
# 获取唯一的标签值
unique_labels = np.unique(labels)
# 创建图例标签
legend_labels = [f"g={int(label)}" for label in unique_labels]
# 创建图例
plt.legend(fontsize=16,bbox_to_anchor=(1.25, 1),loc="upper right",handles=scatter.legend_elements()[0], labels=legend_labels)
plt.subplots_adjust(right=0.8)
plt.title('KMeans Clustering with TSNE',fontsize=25)
plt.savefig(f"tsne.png")

#plt.xlabel('Principal Component 1')
#plt.ylabel('Principal Component 2')
#sys.exit()

# # 创建一个字典来存储每个聚类的索引
# cluster_indices = defaultdict(list)

# for idx, label in enumerate(labels):
#     cluster_indices[label].append(idx)

# 内部评价
#silhouette_avg = silhouette_score(all_rhocurves, labels)
#print("轮廓系数:", silhouette_avg)
print("各类指标:", all_rhocurves_clusters_index)

# 第k类
for g in range(gnumber) :
    random_num = list(random.sample(range(len(all_rhocurves_clusters_index[g] )), 11))
    #random_num = [0]
    random_index = [all_rhocurves_clusters_index[g][j] for j in random_num]
    for i,index in enumerate(random_index) :
        fsize = 30
        #print(index)
        
        plt.figure(figsize=(15, 9)) 
        plt.plot(B_array, all_rhocurves[index,0:len(B_array)],linewidth=3,linestyle='-', color='blue',label="+")
        plt.plot(B_array, all_rhocurves[index,len(B_array):len(B_array)*2],linewidth=3,linestyle='-', color='red',label="-")
        #plt.plot(B_array[0:len(B_array)-1], all_rhocurves_derivative[index,0:len(B_array)-1],linewidth=3,linestyle='-', color='blue',label="+")
        #plt.plot(B_array[0:len(B_array)-1], all_rhocurves_derivative[index,len(B_array)-1:len(B_array)*2-1],linewidth=3,linestyle='-', color='red',label="-")
        plt.xlabel('B(T)',fontsize=fsize)  # x 轴标签
        plt.xticks(fontsize=fsize)
        plt.ylabel(r'$\rho_{xx}$($a.u.$)',fontsize=fsize)  # y 轴标签
        plt.yticks(fontsize=fsize)
        plt.title(f'n1_{curve_param_array[index,0]:.3e}mu1_{curve_param_array[index,1]:.3e}'+\
                f'n2_{curve_param_array[index,2]:.3e}mu2_{curve_param_array[index,3]:.3e}sigmaA_{curve_param_array[index,4]:.1e}', fontsize=20, pad=30)  # 图像标题
        plt.legend(fontsize=18, loc='upper right', bbox_to_anchor=(1.2, 1))  # 显示图例
        plt.subplots_adjust(right=0.85)
        plt.tick_params(pad=8)
        #plt.xlim(0,0.5)
        #plt.xlim(0,0.5)
        plt.ylim(0,1.1)
        # 显示图像
        plt.grid(True)  # 添加网格
        plt.savefig(f"g={g}_{index}.png")
        plt.close()


###------------------- 聚类 -------------------###






