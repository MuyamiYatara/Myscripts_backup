import numpy as np
from joblib import Parallel, delayed
import torch
import modellib
import figplot 


#!---------------  训练单个模型并适当间隔输出结果画图  ---------------!#
def train_single_model(training_model, X_data, target_data, 
                      lr=0.01, max_iter=2000, print_interval=100, jobid=0):
    """
    training_model: 需要训练的模型
    X_data: X轴数据 (numpy数组)
    target_data: 目标曲线的对应Y轴值 (numpy数组)
    lr: 学习率
    max_iter: 最大迭代步数
    print_interval: 间隔多少步输出一次当前结果
    jobid: 用于管理不同参数的输出
    """
    # 转换数据为tensor
    X_tensor = torch.tensor(X_data, dtype=torch.float32)
    target_tensor = torch.tensor(target_data, dtype=torch.float32)
    

    # 初始化模型
    optimizer = torch.optim.Adam(training_model.parameters(), lr=lr)
    
    # 存储训练过程
    loss_history = []
    
    for step in range(max_iter):
        # 前向传播
        pred = training_model(X_tensor)
        loss = torch.mean((pred - target_tensor)**2)
        
        # 反向传播
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
        # 记录损失
        loss_history.append(loss.item())
        
        # 定期输出参数值并画图
        if step % print_interval == 0 or step == max_iter-1:

            params_now = training_model.params.detach().numpy()
            with torch.no_grad():
                Y_now = training_model(X_tensor).numpy()
                
                figplot.compare_target_fit(X_data, target_data, Y_now, figname=f"{jobid}_res_step={step}.png")
            
            
            
            outname = f"{jobid}_train_history"

            # 结果字符串
            result_str = (f"Step {step:4d} | Loss: {loss:.3e} | "
                        f"Params: n1={params_now[0]:.3e} μ1={params_now[1]:.3e} "
                        f"n2={params_now[2]:.3e} μ2={params_now[3]:.3e} "
                        r"\sigma_A" + f"={params_now[4]:.3e}")

            # 写入文件
            with open(outname, "a") as f:
                f.write(result_str + "\n")
            # print(f"Step {step:4d} | Loss: {loss.item():.3e} | "
            #       f"Params: n1={params_now[0]:.3e} μ1={params_now[1]:.3e} "
            #       f"n2={params_now[2]:.3e} μ2={params_now[3]:.3e}"+r"\sigma_A"+f"={params_now[4]:.3e}")
    
    return training_model.params.detach().numpy(), loss_history[-1]


#!---------------  磁阻以及霍尔曲线的n带模型拟合  ---------------!#
def magresistance_fit(modelname, B_array, target_data, initial_params_list,
                         lr=0.01, max_iter=2000, print_interval=100) :
    """
    不返回任何数值,绘制模型训练结束后的参数下得到的曲线与目标曲线的对比, 同时每隔print_interval步会输出一遍当前参数和当前曲线



    modelname: 字符串, 选择输入的是rhoxx还是rhoxy以及用几带模型拟合

        '2b_rhoxx' : modellib.magrhoxx_TwoBandModel,
        '2b_rhoxy' : modellib.magrhoxy_TwoBandModel,
        '3b_rhoxx' : modellib.magrhoxx_ThreeBandModel,
        '3b_rhoxy' : modellib.magrhoxy_ThreeBandModel

    initial_params_list: 列表, 内有多组初始参数组合, 格式如下, 注意每个模型的参数不一样:
    [
        [参数组1],
        [参数组2],
        ...
    ]

    """

    # 生成模拟数据（示例参数）
    #B_array = np.linspace(-10, 10, 201)
    #real_params = [10.0, 0.5, 1.0, -0.3, 0.000]  # 真实参数
    
    # 创建目标数据（需要实际模型计算）
    #with torch.no_grad():
    #    true_model = modellib.magrhoxx_TwoBandModel(real_params)
    #    target_data = true_model(torch.tensor(B_array)).numpy()
    
    # 设置多组初始参数
    #initial_params_list = [
    #    [5.0, -0.5, 2.0, -0.1, 0.00],   # 合理范围初始化
    #]

    n_jobs = len(initial_params_list) #不同初始参数并行训练
    
    modelname_dict = {
        '2b_rhoxx' : modellib.magrhoxx_TwoBandModel,
        '2b_rhoxy' : modellib.magrhoxy_TwoBandModel,
        '3b_rhoxx' : modellib.magrhoxx_ThreeBandModel,
        '3b_rhoxy' : modellib.magrhoxy_ThreeBandModel
    } 

    initial_model_list = [modelname_dict[modelname](i) for i in initial_params_list]

    # 并行训练
    results = Parallel(n_jobs=n_jobs)(
        delayed(train_single_model)(imodel, B_array, target_data, 
        max_iter=max_iter, lr=lr, print_interval=print_interval, jobid=jid)
        for jid,  imodel in enumerate(initial_model_list)
    )
    
    # 选择最佳参数
    best_idx = np.argmin([r[1] for r in results])
    best_params = results[best_idx][0]
    
    print(f"\n最佳参数: n1={best_params[0]:.3e} μ1={best_params[1]:.3e} "
          f"n2={best_params[2]:.3e} μ2={best_params[3]:.3e}"+r"\sigma_A"+f"={best_params[4]:.3e}")
    
    with torch.no_grad():
        best_model = modelname_dict[modelname](best_params)
        fit_data = best_model(torch.tensor(B_array)).numpy()

    figplot.compare_target_fit(B_array, target_data, fit_data, "bestcurve.png")

    return