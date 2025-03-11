import numpy as np
import torch

#!--------------- Global Parameters ---------------!#

electron = 1.602176634e-19  # C



#!--------------- Boltzmann Transport Theory for Resistivity in a Two-Band Model ---------------!#
class magrhoxx_TwoBandModel(torch.nn.Module):
    def __init__(self, init_params):
        """
        init_params: 初始参数张量 [n1, mu1, n2, mu2, sigmaA]
        """
        super().__init__()
        self.params = torch.nn.Parameter(torch.tensor(init_params, dtype=torch.float32))
        

    def forward(self, B):
        """
        B: 磁场强度张量 (形状 [N])
        返回: rho_xx = sigma_inv[0,0] 预测值 (形状 [N])
        """
        
        
        # 解包参数
        n1, mu1, n2, mu2, sigmaA = self.params
        
        # print(list(self.params.detach().numpy()))
        # 参数放缩回正常的单位
        n1 = n1 * 1e16 
        n2 = n2 * 1e16 


        # 计算基础电导率
        sigma1 = electron * n1 * torch.abs(mu1)
        sigma2 = electron * n2 * torch.abs(mu2)
        
        # 为每个B值创建矩阵
        B_sq = B**2
        mu1_sq = mu1**2
        mu2_sq = mu2**2
        
        # 第一个电导矩阵分量（形状 [N, 2, 2]）
        denom1 = 1 + mu1_sq * B_sq
        sigma1_xx = sigma1 / denom1
        sigma1_xy = sigma1 * mu1 * B / denom1
        mat1 = torch.stack([
            torch.stack([sigma1_xx, sigma1_xy], dim=-1),
            torch.stack([-sigma1_xy, sigma1_xx], dim=-1)
        ], dim=-2)

        
        # 第二个电导矩阵分量
        denom2 = 1 + mu2_sq * B_sq
        sigma2_xx = sigma2 / denom2
        sigma2_xy = sigma2 * mu2 * B / denom2
        mat2 = torch.stack([
            torch.stack([sigma2_xx, sigma2_xy], dim=-1),
            torch.stack([-sigma2_xy, sigma2_xx], dim=-1)
        ], dim=-2)
        
        # 反常霍尔分量
        matA = torch.zeros_like(mat1)
        matA[..., 0, 1] = sigmaA*0.0
        matA[..., 1, 0] = -sigmaA*0.0
        
        # 总电导率矩阵
        sigma_total = mat1 + mat2 + matA
        
        # 矩阵求逆并取[0,0]元素, 乘以1e7
        sigma_inv = torch.linalg.inv(sigma_total)
        # print(sigma_inv[..., 0, 0])
        return sigma_inv[..., 0, 0]  # 形状 [N]
    
class magrhoxy_TwoBandModel(torch.nn.Module):
    def __init__(self, init_params):
        """
        init_params: 初始参数张量 [n1, mu1, n2, mu2, sigmaA]
        """
        super().__init__()
        self.params = torch.nn.Parameter(torch.tensor(init_params, dtype=torch.float32))
        
    def forward(self, B):
        """
        B: 磁场强度张量 (形状 [N])
        返回: rho_xy = sigma_inv[0,1] 预测值 (形状 [N])
        """
        # 解包参数
        n1, mu1, n2, mu2, sigmaA = self.params

        # 参数放缩回正常的单位
        n1 = n1 * 1e27 
        n2 = n2 * 1e27 
        
        # 计算基础电导率
        sigma1 = electron * n1 * torch.abs(mu1)
        sigma2 = electron * n2 * torch.abs(mu2)
        
        # 为每个B值创建矩阵
        B_sq = B**2
        mu1_sq = mu1**2
        mu2_sq = mu2**2
        
        # 第一个电导矩阵分量（形状 [N, 2, 2]）
        denom1 = 1 + mu1_sq * B_sq
        sigma1_xx = sigma1 / denom1
        sigma1_xy = sigma1 * mu1 * B / denom1
        mat1 = torch.stack([
            torch.stack([sigma1_xx, sigma1_xy], dim=-1),
            torch.stack([-sigma1_xy, sigma1_xx], dim=-1)
        ], dim=-2)
        
        # 第二个电导矩阵分量
        denom2 = 1 + mu2_sq * B_sq
        sigma2_xx = sigma2 / denom2
        sigma2_xy = sigma2 * mu2 * B / denom2
        mat2 = torch.stack([
            torch.stack([sigma2_xx, sigma2_xy], dim=-1),
            torch.stack([-sigma2_xy, sigma2_xx], dim=-1)
        ], dim=-2)
        
        # 反常霍尔分量
        matA = torch.zeros_like(mat1)
        matA[..., 0, 1] = sigmaA
        matA[..., 1, 0] = -sigmaA
        
        # 总电导率矩阵
        sigma_total = mat1 + mat2 + matA
        
        # 矩阵求逆并取[0,0]元素
        sigma_inv = torch.linalg.inv(sigma_total)  # 形状 [N, 2, 2]
        return sigma_inv[..., 0, 1]   # 形状 [N]


class magrhoxx_ThreeBandModel(torch.nn.Module):
    def __init__(self, init_params):
        """
        init_params: 初始参数张量 [n1, mu1, n2, mu2, n3, mu3, sigmaA]
        """
        super().__init__()
        self.params = torch.nn.Parameter(torch.tensor(init_params, dtype=torch.float32))
        

    def forward(self, B):
        """
        B: 磁场强度张量 (形状 [N])
        返回: rho_xx = sigma_inv[0,0] 预测值 (形状 [N])
        """
        
        
        # 解包参数
        n1, mu1, n2, mu2, n3, mu3, sigmaA = self.params
        
        # print(list(self.params.detach().numpy()))
        # 参数放缩回正常的单位
        n1 = n1 * 1e15 
        n2 = n2 * 1e15
        n3 = n3 * 1e15 


        # 计算基础电导率
        sigma1 = electron * n1 * torch.abs(mu1)
        sigma2 = electron * n2 * torch.abs(mu2)
        sigma3 = electron * n3 * torch.abs(mu3)
        
        # 为每个B值创建矩阵
        B_sq = B**2
        mu1_sq = mu1**2
        mu2_sq = mu2**2
        mu3_sq = mu3**2
        
        # 第一个电导矩阵分量（形状 [N, 2, 2]）
        denom1 = 1 + mu1_sq * B_sq
        sigma1_xx = sigma1 / denom1
        sigma1_xy = sigma1 * mu1 * B / denom1
        mat1 = torch.stack([
            torch.stack([sigma1_xx, sigma1_xy], dim=-1),
            torch.stack([-sigma1_xy, sigma1_xx], dim=-1)
        ], dim=-2)

        
        # 第二个电导矩阵分量
        denom2 = 1 + mu2_sq * B_sq
        sigma2_xx = sigma2 / denom2
        sigma2_xy = sigma2 * mu2 * B / denom2
        mat2 = torch.stack([
            torch.stack([sigma2_xx, sigma2_xy], dim=-1),
            torch.stack([-sigma2_xy, sigma2_xx], dim=-1)
        ], dim=-2)

        # 第三个电导矩阵分量
        denom3 = 1 + mu3_sq * B_sq
        sigma3_xx = sigma3 / denom3
        sigma3_xy = sigma3 * mu3 * B / denom3
        mat3 = torch.stack([
            torch.stack([sigma3_xx, sigma3_xy], dim=-1),
            torch.stack([-sigma3_xy, sigma3_xx], dim=-1)
        ], dim=-2)
        
        # 反常霍尔分量
        matA = torch.zeros_like(mat1)
        matA[..., 0, 1] = sigmaA*0.0
        matA[..., 1, 0] = -sigmaA*0.0
        
        # 总电导率矩阵
        sigma_total = mat1 + mat2 + mat3 + matA
        
        # 矩阵求逆并取[0,0]元素, 乘以1e7
        sigma_inv = torch.linalg.inv(sigma_total)
        # print(sigma_inv[..., 0, 0])
        return sigma_inv[..., 0, 0]  # 形状 [N]

class magrhoxy_ThreeBandModel(torch.nn.Module):
    def __init__(self, init_params):
        """
        init_params: 初始参数张量 [n1, mu1, n2, mu2, n3, mu3, sigmaA]
        """
        super().__init__()
        self.params = torch.nn.Parameter(torch.tensor(init_params, dtype=torch.float32))
        

    def forward(self, B):
        """
        B: 磁场强度张量 (形状 [N])
        返回: rho_xx = sigma_inv[0,0] 预测值 (形状 [N])
        """
        
        
        # 解包参数
        n1, mu1, n2, mu2, n3, mu3, sigmaA = self.params
        
        # print(list(self.params.detach().numpy()))
        # 参数放缩回正常的单位
        n1 = n1 * 1e16 
        n2 = n2 * 1e16
        n3 = n3 * 1e16 


        # 计算基础电导率
        sigma1 = electron * n1 * torch.abs(mu1)
        sigma2 = electron * n2 * torch.abs(mu2)
        sigma3 = electron * n3 * torch.abs(mu3)
        
        # 为每个B值创建矩阵
        B_sq = B**2
        mu1_sq = mu1**2
        mu2_sq = mu2**2
        mu3_sq = mu3**2
        
        # 第一个电导矩阵分量（形状 [N, 2, 2]）
        denom1 = 1 + mu1_sq * B_sq
        sigma1_xx = sigma1 / denom1
        sigma1_xy = sigma1 * mu1 * B / denom1
        mat1 = torch.stack([
            torch.stack([sigma1_xx, sigma1_xy], dim=-1),
            torch.stack([-sigma1_xy, sigma1_xx], dim=-1)
        ], dim=-2)

        
        # 第二个电导矩阵分量
        denom2 = 1 + mu2_sq * B_sq
        sigma2_xx = sigma2 / denom2
        sigma2_xy = sigma2 * mu2 * B / denom2
        mat2 = torch.stack([
            torch.stack([sigma2_xx, sigma2_xy], dim=-1),
            torch.stack([-sigma2_xy, sigma2_xx], dim=-1)
        ], dim=-2)

        # 第三个电导矩阵分量
        denom3 = 1 + mu3_sq * B_sq
        sigma3_xx = sigma3 / denom3
        sigma3_xy = sigma3 * mu3 * B / denom3
        mat3 = torch.stack([
            torch.stack([sigma3_xx, sigma3_xy], dim=-1),
            torch.stack([-sigma3_xy, sigma3_xx], dim=-1)
        ], dim=-2)
        
        # 反常霍尔分量
        matA = torch.zeros_like(mat1)
        matA[..., 0, 1] = sigmaA
        matA[..., 1, 0] = -sigmaA
        
        # 总电导率矩阵
        sigma_total = mat1 + mat2 + mat3 + matA
        
        # 矩阵求逆并取[0,0]元素, 乘以1e7
        sigma_inv = torch.linalg.inv(sigma_total)
        # print(sigma_inv[..., 0, 0])
        return sigma_inv[..., 0, 1]  # 形状 [N]

if __name__ == "__main__":
    
    B_array = np.linspace(-5, 5, 100)
    real_params = [1e28, 0.5, 1e27, -0.3, 1e3]
    with torch.no_grad():
        rmat1 = magrhoxx_TwoBandModel(real_params)
        print(rmat1(torch.tensor(B_array, dtype=torch.float32)).numpy())

