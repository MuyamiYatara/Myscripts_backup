import numpy as np
import itertools
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from matplotlib.colors import LinearSegmentedColormap
#a = np.array([1,2,3,4,5])
#b = np.array([1,2,3,4])
#c = np.array([1,2,3])
#d = np.array([1,2])
#
#grid1, grid2, grid3, grid4 = np.meshgrid(a,b,c,d)
#result = np.vstack([grid1.ravel(), grid2.ravel(), grid3.ravel(), grid4.ravel()]).T
#print(np.array(list(itertools.product(a,b,c,d))))

kb = 1.380649*10**-23
mub = 9.28*10**-24
T = 100.0
Tc = 175.0 
def equation(x):
    return x - np.tanh( mub*1.0/(kb*T)* B_v + (Tc / T) * x)


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

fsize = 30
B_array = [-1 + j*(2)/200 for j in range(200)]
T_array = [20 + j*(160)/10 for j in range(10)]

plt.figure(figsize=(15, 9)) 

num_colors = 10
colors = ["red", "blue"]
cmap_name = "o_blue"
blue_red_cmap = LinearSegmentedColormap.from_list(cmap_name, colors)
n = 10  # 设置曲线数量
colors = [blue_red_cmap(i / (n )) for i in range(n)]


# for ind, i in enumerate(B_array) :
#     Mcuv = []
#     B_v = i
#     for t in T_array :
#         T = t
    
#         sols = fsolve(equation, 0.7)
#         Mcuv.append(sols[0]*np.tanh(5.0*B_v))
#         #Mcuv.append(np.tanh(0.8*B_v))
#     plt.plot(T_array, Mcuv,linewidth=3,linestyle='-', color=colors[ind],label=f"B={i}T")
# plt.xlabel('T(K)',fontsize=fsize)  # x 轴标签




# for ind, t in enumerate(T_array) :
#     Mcuv = []
#     T = t
#     for i in B_array :
#         B_v = i
#         if i>0 :
#             sols = fsolve(equation, 0.7)
#         else :
#             sols = -fsolve(equation, -0.7)
#         Mcuv.append(sols[0]*np.tanh(5*B_v))
#         #Mcuv.append(np.tanh(0.8*B_v))
#     if ind == 0 :
#         #for l in range(15, len(B_array)) :
#         l = 10
#         #Mcuv_shifted_minus = np.array(Mcuv)-(Mcuv[l]-Mcuv[0])
#         B_array_shifted_minus = np.array(B_array)-(B_array[l]-B_array[0])
#         #Mcuv_shifted_posit = np.array(Mcuv)+(Mcuv[l]-Mcuv[0])
#         B_array_shifted_posit = np.array(B_array)+(B_array[l]-B_array[0])

#         dk = (Mcuv[l+1] - Mcuv[l])/(B_array[l+1] - B_array[l])
#         b1 = Mcuv[l] - B_array[l]*dk
#         b2 = Mcuv[len(Mcuv)-l-1] - B_array[len(B_array)-l-1]*dk
#         left_y =  dk*B_array[0] + b1
#         right_y = dk*B_array[len(B_array)-1] + b2
         
#         #print(l, np.abs(B_array_shifted_posit[l]*dk+b1 - Mcuv_shifted_posit[l]))
#         #if np.abs(B_array_shifted_posit[l]*dk+b1 - Mcuv_shifted_posit[l]) <= 0.1 :
#         #    break

#         #plt.plot(B_array, np.array(B_array)*dk+b1,linewidth=3,linestyle='-', color="black",label=f"line1")
#         #plt.plot(B_array, np.array(B_array)*dk+b2,linewidth=3,linestyle='-', color="green",label=f"line2")
#         #plt.plot(np.array(B_array)-(B_array[l]-B_array[0]), np.array(Mcuv)-(Mcuv[l]-left_y),linewidth=3,linestyle='-', color="red",label=f"-1")
#         #plt.plot(np.array(B_array)-(B_array[len(B_array)-1-l]-B_array[len(B_array)-1]), np.array(Mcuv)-(Mcuv[len(B_array)-1-l]-right_y),linewidth=3,linestyle='-', color="blue",label=f"+1")

        
#         plt.plot(np.array(B_array[0:len(B_array)-l]) - (B_array[0] - B_array[l]), np.array(Mcuv[0:len(B_array)-l])-(Mcuv[l]-left_y),linewidth=3,linestyle='-', color="red",label=f"-1")
#         plt.plot(np.array(B_array[l:len(B_array)]) - (B_array[len(B_array)-1] - B_array[len(B_array)-l-1]), np.array(Mcuv[l:len(B_array)])-(Mcuv[len(B_array)-1-l]-right_y),linewidth=3,linestyle='-', color="blue",label=f"+1")
#         plt.plot(B_array[0:l], np.array(B_array[0:l])*dk+b1,linewidth=3,linestyle='-', color="black",label=f"line1")
#         plt.plot(B_array[len(B_array)-l:len(B_array)], np.array(B_array[len(B_array)-l:len(B_array)])*dk+b2,linewidth=3,linestyle='-', color="green",label=f"line2")



loop1 = M_BT_curve(B_array, 100.0, 5.0, 100)

plt.plot(B_array, loop1[0:len(B_array)],linewidth=3,linestyle='-', color="red",label=f"-1")
plt.plot(B_array, loop1[len(B_array):len(B_array)*2],linewidth=3,linestyle='-', color="red",label=f"-1")
plt.xlabel('B(T)',fontsize=fsize)  # x 轴标签



plt.xticks(fontsize=fsize)
plt.ylabel(r'$M$/$M_s$',fontsize=fsize)  # y 轴标签
plt.yticks(fontsize=fsize)
plt.legend(fontsize=18, loc='upper right', bbox_to_anchor=(1.2, 1))  # 显示图例
plt.subplots_adjust(right=0.85)
plt.tick_params(pad=8)
#plt.xlim(0,0.5)
#plt.ylim(0,40)
# 显示图像
plt.grid(True)  # 添加网格
plt.savefig(f"M.png")


