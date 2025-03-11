
import torch
import numpy as np
import gradfit 
import figplot, readfile, modellib


# 5. 主程序
if __name__ == "__main__":

    # mat_name = "../3V_150535_ScanB_I7-1Vxx5-3&8-10Vxy10-3&8-5ground2-11_hengya100-10k_RI1k__T1.9992K_Vbg3V_Vbias0V_ac1V.mat"
    # B_and_target = readfile.read_matlab(mat_name, 'B', 'Rxx')
    # init_p = [
    #     [1.3, 0.16, 1.1, -0.05, 0.0],
    # ]
    mat_name = "../3V_150535_ScanB_I7-1Vxx5-3&8-10Vxy10-3&8-5ground2-11_hengya100-10k_RI1k__T1.9992K_Vbg3V_Vbias0V_ac1V.mat"
    B_and_target = readfile.read_matlab(mat_name, 'B', 'Rxx')
    init_p = [
        np.array([1.08, 0.25, 4.51, -1.5, 0.59, 2.0, 0.0])/2,
    ]
    gradfit.magresistance_fit('3b_rhoxx', B_and_target['B']+0.001, B_and_target['Rxx'], initial_params_list=init_p ,lr=0.01, max_iter=10000, print_interval=1000)
