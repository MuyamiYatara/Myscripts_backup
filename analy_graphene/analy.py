import numpy as np
import matplotlib.pyplot as plt
from scipy.special import digamma


#------------------------ 常量定义 ------------------------#
# 物理常数
kb = 1.380649*10**-23    # J/K
e = 1.602177*10**-19     # C
vf = 0.8*10**6               # m/s
hbar = 1.054572*10**-34                  # J*s
#化学势 (eV) 
#mu = -0.0

#tau公式内的系数
tau0 = 2
alpha = 8*10**(-5)
beta = 0*2*10**(-13)
thetaD = 300.0
A = 8*10**3

#设置字体
title_font_size = 23
tick_font_size = 21
lw = 4


#------------------------ 常量定义 ------------------------#









#------------------------ 函数定义 ------------------------#
def dfermi(E, mu, T) :
    return 1.0/(kb*T/e) * (np.exp((E-mu)/(kb*T/e)))/(np.exp((E-mu)/(kb*T/e)) + 1.0)**2

def intgrand_ne(E, mu, T) :
    return 1/(np.exp((E-mu)/(kb*T/e)) + 1.0)*(E/(vf**2))*(e/hbar)**2*(2/np.pi)

def intgrand_nh(E, mu, T) :
    return -np.exp((E-mu)/(kb*T/e))/(np.exp((E-mu)/(kb*T/e)) + 1.0)*(E/(vf**2))*(e/hbar)**2*(2/np.pi)

def intgrand_xx(E, mu, T, B, tau):
    #tau = (1/(tau0+T**2*alpha + T**5*beta))
    Btau = B*tau
    return dfermi(E, mu, T)*(E**2*np.abs(E))/(E**2 + (e*vf**2*Btau/e)**2)*(e/hbar)**2*(1/e)*e**2*(1/(2*np.pi))

def intgrand_xy(E, mu, T, B, tau):
    #tau = (1/(tau0+T**2*alpha + T**5*beta))
    Btau = B*tau
    return -dfermi(E, mu, T)*(E*np.abs(E)*(e*Btau*vf**2/e))/(E**2 + (e*vf**2*Btau/e)**2)*(e/hbar)**2*(1/e)*e**2*(1/(2*np.pi))

def rho_xx_calc(sigma_xx, sigma_xy):
    return sigma_xx/(sigma_xx**2 + sigma_xy**2)

def tau_intgrand(y):
    return (y**5/((np.exp(y) - 1)*(1 - np.exp(-y))))

def F_func(z):
    return np.log(z) + digamma(0.5 + 1.0/z)

def sigma_xx_McCann(BTau, T):
    tau_phi = 10.0/(1+T/10)
    tau_e = 5
    tau_star = 1.0
    return (F_func(2*e*vf**2*BTau/tau_phi) - F_func(2*e*vf**2*BTau/(tau_phi + tau_e)) - 2*F_func(2*e*vf**2*BTau/(tau_phi + tau_star)))*(e**2)


#------------------------ 函数定义 ------------------------#








#------------------------ 计算 ------------------------#
# 积分区间
a, b = -0.5, 0.50
x = np.linspace(a, b, 5000)  # 取足够多的点
x_ne = np.linspace(0.0, b, 5000)
x_nh = np.linspace(a, 0.0, 5000)
#生成btau与T数组
bt_count = 10
t_count = 300
mu_count = 6
B_array = np.linspace(0, 9, bt_count)
Btau_array = np.linspace(0, 9*10**-13, bt_count)
T_array = np.linspace(20, 300, t_count)
mu_array = np.linspace(0.0, 0.0001, mu_count)


ne_array = np.zeros((t_count, mu_count))
nh_array = np.zeros((t_count, mu_count))
tau_array = np.zeros((t_count))
rho_xx_mu_array = np.zeros((bt_count, t_count, mu_count, ))
for m, mu in enumerate(mu_array) :



    sigma_xx_array = np.zeros((bt_count, t_count))
    #sigma_xx_with_WL_array = np.zeros((bt_count, t_count))
    sigma_xy_array = np.zeros((bt_count, t_count))

    sigma_xx_tau_array = np.zeros((bt_count, t_count))
    sigma_xy_tau_array = np.zeros((bt_count, t_count))

    rho_xx_array = np.zeros((bt_count, t_count))
    rho_xx_tau_array = np.zeros((bt_count, t_count))
    #rho_xx_with_WL_array = np.zeros((bt_count, t_count))
    # rho_B0_array = np.zeros((t_count))
    

    for j, t in enumerate(T_array) :
        y = np.linspace(0.0001, thetaD/t, 5000)
        tau_T = tau_intgrand(y)
        #print(tau_T)
        rho_B0 = (np.trapz(tau_T, y)*A*(t/thetaD)**5 + 2*10**3)
        
        sigma_xx_B0 = np.trapz(intgrand_xx(x, 0.002, t, 0.0, 0.0), x)
        sigma_xy_B0 = np.trapz(intgrand_xy(x, 0.002, t, 0.0, 0.0), x)
        rhotau_B0 = rho_xx_calc(sigma_xx_B0, sigma_xy_B0)
        tau = rhotau_B0/rho_B0
        print(tau, rhotau_B0, rho_B0)
        tau_array[j] = tau
        
        ne_array[j, m] = np.trapz(intgrand_ne(x_ne, mu, t), x_ne)
        nh_array[j, m] = np.trapz(intgrand_nh(x_nh, mu, t), x_nh)



        for i, bt in enumerate(B_array) :
        
            #sxx_MC = sigma_xx_McCann(bt, t)
            s_xx = intgrand_xx(x, mu, t, bt, tau)
            s_xy = intgrand_xy(x, mu, t, bt, tau)
            s_tau_xx = intgrand_xx(x, mu, t, Btau_array[i], tau=1.0)
            s_tau_xy = intgrand_xy(x, mu, t, Btau_array[i], tau=1.0)
            

            sigma_xx_tau_array[i, j] = np.trapz(s_tau_xx, x)
            sigma_xy_tau_array[i, j] = np.trapz(s_tau_xy, x)
            sigma_xx_array[i, j] = np.trapz(s_xx, x)*tau
            #sigma_xx_with_WL_array[i, j] = sxx_MC
            sigma_xy_array[i, j] = np.trapz(s_xy, x)*tau


            rho_xx_array[i, j] = rho_xx_calc(sigma_xx_array[i, j], sigma_xy_array[i, j])
            rho_xx_mu_array[i, j, m] = rho_xx_array[i, j]
            rho_xx_tau_array[i, j] = rho_xx_calc(sigma_xx_tau_array[i, j], sigma_xy_tau_array[i, j])
            #rho_xx_with_WL_array[i, j] = rho_xx_calc(sigma_xx_with_WL_array[i, j], sigma_xy_array[i, j])



    #------------------------ 计算 ------------------------#






    #------------------------ 绘图 ------------------------#


    # # df(E)/dE-E
    # plt.figure(figsize=(10,8))
    # for j, T in enumerate(T_array):
    #     plt.plot(x, dfermi(x, mu, T), label=f"T={T:.0f}", linewidth=lw)

    # plt.tick_params(axis='both', labelsize=tick_font_size)
    # plt.xlabel("E",fontsize=tick_font_size)
    # plt.ylabel(r"df(E)/dE",fontsize=tick_font_size)
    # plt.title(r"df(E)/dE vs E for different T",fontsize=title_font_size)
    # plt.legend(fontsize=tick_font_size)
    # plt.grid(True)
    # plt.tight_layout()
    # plt.savefig("dfermi.png")





    # # sigma_xx-Btau
    # plt.figure(figsize=(10,8))

    # for j, T in enumerate(T_array):
    #     plt.plot(B_array, sigma_xx_array[:, j], label=f"noWL T={T:.0f}", linewidth=lw)
    #     plt.plot(B_array, sigma_xx_with_WL_array[:, j], label=f"WL T={T:.0f}", linewidth=lw)

    # # for j, T in enumerate(B_array):
    # #     plt.plot(T_array, sigma_xx_array[i, :], label=f"btau={T:.0f}")

    # #for j, T in enumerate(T_array):
    # #    plt.plot(x, intgrand_xx(x, mu, T, bt), label=f"T={T:.0f}")
    # plt.tick_params(axis='both', labelsize=tick_font_size)
    # plt.xlabel("BTau",fontsize=tick_font_size)
    # plt.ylabel(r"$\sigma_{xx}$",fontsize=tick_font_size)
    # plt.title(r"$\sigma_{xx}$ vs bTau for different T",fontsize=title_font_size)
    # plt.legend(fontsize=tick_font_size)
    # plt.grid(True)
    # plt.tight_layout()
    # plt.savefig("sigma_xx.png")



    # # sigma_xy-Btau
    # plt.figure(figsize=(10,8))

    # for j, T in enumerate(T_array):
    #     plt.plot(B_array, sigma_xy_array[:, j], label=f"T={T:.0f}", linewidth=lw)

    # # for j, T in enumerate(B_array):
    # #     plt.plot(T_array, sigma_xx_array[i, :], label=f"btau={T:.0f}")

    # #for j, T in enumerate(T_array):
    # #    plt.plot(x, intgrand_xx(x, mu, T, bt), label=f"T={T:.0f}")

    # plt.tick_params(axis='both', labelsize=tick_font_size)
    # plt.xlabel("BTau",fontsize=tick_font_size)
    # plt.ylabel(r"$\sigma_{xy}$",fontsize=tick_font_size)
    # plt.title(r"$\sigma_{xy}$ vs bTau for different T",fontsize=title_font_size)
    # plt.legend(fontsize=tick_font_size)
    # plt.grid(True)
    # plt.tight_layout()
    # plt.savefig("sigma_xy.png")



  


    # # rho_xx-T
    # plt.figure(figsize=(10,8))

    # for j, bt in enumerate(B_array):
    #     plt.plot(T_array, rho_xx_array[j, :], label=f"B={bt:.0f}", linewidth=lw)
        

    # # for j, T in enumerate(B_array):
    # #     plt.plot(T_array, sigma_xx_array[i, :], label=f"btau={T:.0f}")

    # #for j, T in enumerate(T_array):
    # #    plt.plot(x, intgrand_xx(x, mu, T, bt), label=f"T={T:.0f}")
    # #plt.ylim(0, 10)
    # plt.tick_params(axis='both', labelsize=tick_font_size)
    # plt.xlabel("T/K",fontsize=tick_font_size)
    # plt.ylabel(r"$\rho_{xx}$",fontsize=tick_font_size)
    # plt.title(r"$\rho_{xx}$ vs T for different B with $\mu$ = "+f"{mu:.4f}",fontsize=title_font_size)
    # plt.legend(fontsize=tick_font_size)
    # plt.grid(True)
    # plt.tight_layout()
    # plt.savefig(f"rho_T_mu_{mu}.png")

# rho_xx-T
plt.figure(figsize=(10,8))

for m, mu in enumerate(mu_array):
    if (m == 0) : continue
    plt.plot(T_array, rho_xx_mu_array[9, :, m], label=f"mu={mu:.6f}", linewidth=lw)
        

    # for j, T in enumerate(B_array):
    #     plt.plot(T_array, sigma_xx_array[i, :], label=f"btau={T:.0f}")

    #for j, T in enumerate(T_array):
    #    plt.plot(x, intgrand_xx(x, mu, T, bt), label=f"T={T:.0f}")
    #plt.ylim(0, 10)
plt.tick_params(axis='both', labelsize=tick_font_size)
plt.xlabel("T/K",fontsize=tick_font_size)
plt.ylabel(r"$\rho_{xx}$",fontsize=tick_font_size)
plt.title(r"$\rho_{xx}$ vs T for different $\mu$ with B  = "+f"9T",fontsize=title_font_size)
plt.legend(fontsize=tick_font_size)
plt.grid(True)
plt.tight_layout()
plt.savefig(f"rho_T-diffmu.png")

# # tau-T
# plt.figure(figsize=(10,8))

# plt.plot(T_array, tau_array[:], linewidth=lw)

# plt.tick_params(axis='both', labelsize=tick_font_size)
# plt.xlabel("T/K",fontsize=tick_font_size)
# plt.ylabel(r"$\tau$",fontsize=tick_font_size)
# plt.title(r"$\tau$ vs T",fontsize=title_font_size)
# #plt.legend(fontsize=tick_font_size)
# plt.grid(True)
# plt.tight_layout()
# plt.savefig(f"tau_T.png")


# # ne-T
# plt.figure(figsize=(10,8))
# for m, mu in enumerate(mu_array):
#     plt.plot(T_array, ne_array[:,m],label=f"mu={mu:.3f}", linewidth=lw)
# plt.ylim(0,2*10**15)
# plt.tick_params(axis='both', labelsize=tick_font_size)
# plt.xlabel("T/K",fontsize=tick_font_size)
# plt.ylabel(r"$n_{e}$/$m^{-2}$",fontsize=tick_font_size)
# plt.title(r"$n_{e}$ vs T",fontsize=title_font_size)
# plt.legend(fontsize=tick_font_size)
# plt.grid(True)
# plt.tight_layout()
# plt.savefig(f"ne_T.png")

# # nh-T
# plt.figure(figsize=(10,8))
# for m, mu in enumerate(mu_array):
#     plt.plot(T_array, nh_array[:,m],label=f"mu={mu:.3f}", linewidth=lw)
# plt.ylim(0,2*10**15)
# plt.tick_params(axis='both', labelsize=tick_font_size)
# plt.xlabel("T/K",fontsize=tick_font_size)
# plt.ylabel(r"$n_{h}$/$m^{-2}$",fontsize=tick_font_size)
# plt.title(r"$n_{h}$ vs T",fontsize=title_font_size)
# plt.legend(fontsize=tick_font_size)
# plt.grid(True)
# plt.tight_layout()
# plt.savefig(f"nh_T.png")



#------------------------ 绘图 ------------------------#