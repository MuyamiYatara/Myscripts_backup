from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.colors import LinearSegmentedColormap



theta = 10.0
rho0 = 0.733
alpha = 50.0

# 定义被积函数
def integrand(x):
    return x**2/((np.exp(x) - 1.0)*(1.0 - np.exp(-x)))

Tlist = np.linspace(30, 300, 31)
#rho_nomag = []



##plt.figure(figsize=(8, 6)) 
##plt.plot(Tlist, rho_nomag, linestyle='-', color='blue', label='bg_formula')
##
### 添加图例和标签
##plt.xlabel('x')  # x 轴标签
##plt.ylabel('y')  # y 轴标签
##plt.title('Example Plot')  # 图像标题
##plt.legend()  # 显示图例
##plt.xlim(0,300)
### 显示图像
##plt.grid(True)  # 添加网格
##plt.savefig("rho_fit.png")



import numpy as np


def read_data(file_name):
    data = {}
    iband = None
    T = None
    block_data = []
    
    # 指定编码方式或忽略解码错误
    with open(file_name, 'r', encoding='utf-8', errors='ignore') as file:
        # 跳过前五行并直接从第六行开始处理
        for _ in range(4):
            next(file)
        
        T = None
        data = {}
        flag = 1
        # 从第六行开始读取
        for line in file:
            line = line.strip()

            #if line.startswith("#  iband"):
            #    if T is not None :
            #        data[iband][T] = np.array(block_data)
            #    # 新的iband块
            #    iband = int(line.split('=')[1].strip())
            #    print(iband)
            #    if iband not in data:
            #        data[iband] = {}
            #        flag = 1
            #        T = None
            #   
            if line.startswith("#  T ="):
                # 存储当前T的数据
                if (flag == 1):
                    flag = 0
                    T = line.split('=')[1].strip()
                    T = float(T.split(' ')[0].strip()) 
                    #print(T)

                if (flag == 0):    
                    data[T] = np.array(block_data)
                    T = line.split('=')[1].strip()
                    T = float(T.split(' ')[0].strip())
                    
                block_data = []
            elif line and not line.startswith("#"):
                # 读取数据行
                row = list(map(float, line.split()))
                block_data.append(row)
                
        # 最后一个T块需要保存
        if T is not None:
            data[T] = np.array(block_data)

    return data

def cal_n_iband(data, ne_filename) :
    with open(ne_filename, 'w') as f:
        for iband in data :
            f.write(f"#  iband =  {iband}\n")
            f.write(f"# Column    1               2\n")
            f.write(f"#           T/K              n/m^-3\n")
            for T in data[iband] :

                abs_sigmaxy = np.abs(data[iband][T][:,2])
                max_row = np.argmax(abs_sigmaxy)
                Btau_max_sigmaxy = np.abs(data[iband][T][max_row,0])
                sigmaxx_0 = data[iband][T][0,1]

                ne_iband = 1.0/(1.602*10**-19)*Btau_max_sigmaxy*sigmaxx_0
                f.write(f"{T}    {ne_iband}\n")
            f.write(f"\n")



# 用法示例
numbers = np.linspace(-0.01000, 0.01000, 11)
formatted_numbers = [f"{num:.5f}" for num in numbers]
mu = '0.00000'
sigma_file_name = './magR/rho_total_mu_'+ mu +'eV.dat'  # 输入文件名
data = read_data(sigma_file_name)
#print(Tlist, "Tlist")
#line111 = '0.000000E+00    0.671615E-18   -0.698571E-34   -0.744137E-34   -0.712976E-34    0.671615E-18   -0.170473E-33   -0.744137E-34   -0.170473E-33    0.131549E-17'
#print(np.array(list(map(float, line111.split()))))
#print(data[30.0])
#ne_file_name = 'ne_bands_mu_'+ n +'eV.dat'
#cal_n_iband(data, ne_file_name)
tau = 0.0
#Tlist = []
Blist = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
Brho = {}
fl = 1
for T,line in data.items() :
    if (fl == 1) :
        fl = 0
        result = 0.0
        # 调用 quad 函数计算积分
        result, error = quad(integrand, 0.0, theta/(T+0.1))
        rho_nomag = rho0 + alpha*(T/theta)**2*result # rho_nomag in  unit of muOhm*cm
        tau =  line[0,1]/(rho_nomag*10**(-8))*10**12 #tau in unit of ps
        #Tlist.append(T)
        # 线性插值
        linear_interp = interp1d(line[:,0]/tau, line[:,1]/(tau*10**-12), kind='linear')
        # 插值点
        #x_new = np.linspace(x[0], x[-1], 100)
        BTval = [linear_interp(b) for b in Blist]
        Brho = np.array(BTval)
        #Brho[line[iB,0]/tau] = line[iB,:]/tau*10**12
        #Brho[line[iB,0]/tau][0] = T         
    else :
        result = 0.0
        #print(111)
        # 调用 quad 函数计算积分
        result, error = quad(integrand, 0.0, theta/(T+0.1))
        rho_nomag = rho0 + alpha*(T/theta)**2*result # rho_nomag in  unit of muOhm*cm
        tau =  line[0,1]/(rho_nomag*10**(-8))*10**12 #tau in unit of ps
        #Tlist.append(T)
        linear_interp = interp1d(line[:,0]/tau, line[:,1]/(tau*10**-12), kind='linear')
        # 插值点
        #x_new = np.linspace(x[0], x[-1], 100)
        BTval = [linear_interp(b) for b in Blist]
        Brho = np.vstack((Brho, np.array(BTval)))

#print(Brho)

colors = ["blue", "red"]
cmap_name = "blue_red"
blue_red_cmap = LinearSegmentedColormap.from_list(cmap_name, colors)
n = 11  # 设置曲线数量
colors = [blue_red_cmap(i / (n - 1)) for i in range(n)]

plt.figure(figsize=(12, 8)) 
ax = plt.gca()
ax.yaxis.get_offset_text().set_fontsize(16)
for i in range(n) :
    plt.plot(Tlist,Brho[:,i] , linestyle='-', color=colors[i], label=f'B={i}')


# 添加图例和标签
plt.xlabel('T(K)',fontsize=19)  # x 轴标签
plt.xticks(fontsize=19)
plt.ylabel(r'$\rho_{xx}$($\Omega \cdot s$)',fontsize=19)  # y 轴标签
plt.yticks(fontsize=19)
plt.title(r'$\rho_{xx}-T$ curve of different B in $\mu =$' + mu+r'$eV$', fontsize=22, pad=30)  # 图像标题
plt.legend(fontsize=18, loc='upper right', bbox_to_anchor=(1.25, 1))  # 显示图例
plt.subplots_adjust(right=0.80)

plt.xlim(0,300)
# 显示图像
plt.grid(True)  # 添加网格
plt.savefig("rho_fit.png")


    



