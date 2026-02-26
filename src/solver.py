# Core algorithm
import torch
import numpy as np

class scSCF_Solver:
    def __init__(self, X_vac, X_tissue, gene_names):
        """
        X_vac: numpy array (n_cells_vac, n_genes)
        X_tissue: numpy array (n_cells_tissue, n_genes)
        """
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        
        # 1. 计算基态统计量 (Vacuum Statistics)
        # 转换为Tensor
        vac_tensor = torch.tensor(X_vac, dtype=torch.float32).to(self.device)
        self.mu_vac = torch.mean(vac_tensor, dim=0)
        
        # 计算协方差矩阵 (Susceptibility Matrix Chi)
        # 注意：对于高维数据，这里可能需要使用PCA降维或Shrinkage Covariance estimator
        vac_centered = vac_tensor - self.mu_vac
        self.C_0 = torch.mm(vac_centered.T, vac_centered) / (vac_tensor.shape[0] - 1)
        
        # 2. 计算观测到的漂移 (Observed Shift)
        tissue_tensor = torch.tensor(X_tissue, dtype=torch.float32).to(self.device)
        self.mu_tissue = torch.mean(tissue_tensor, dim=0)
        self.delta_mu_obs = self.mu_tissue - self.mu_vac
        
        # 3. 初始化待解的场 (Field Parameter)
        # h_env 对应每个基因受到的外力，初始化为0
        self.h_env = torch.zeros(len(gene_names), requires_grad=True, device=self.device)
        
    def define_mask(self, receptors_indices):
        """
        限制场只能作用于受体基因 (Biologically constrained)
        """
        self.mask = torch.zeros_like(self.h_env)
        self.mask[receptors_indices] = 1.0

    def train(self, epochs=1000, lr=0.01, l1_lambda=0.1):
        optimizer = torch.optim.Adam([self.h_env], lr=lr)
        
        for i in range(epochs):
            optimizer.zero_grad()
            
            # 应用掩码：只允许特定基因产生场效应
            effective_h = self.h_env * self.mask if hasattr(self, 'mask') else self.h_env
            
            # 预测漂移 (Linear Response Prediction)
            # delta_mu_pred = C_0 * h
            delta_mu_pred = torch.mv(self.C_0, effective_h)
            
            # 计算损失
            # 重构误差
            mse_loss = torch.norm(self.delta_mu_obs - delta_mu_pred) ** 2
            # 稀疏约束 (L1 Regularization) - 寻找最关键的互作
            sparsity_loss = l1_lambda * torch.norm(effective_h, 1)
            
            loss = mse_loss + sparsity_loss
            
            loss.backward()
            optimizer.step()
            
            if i % 100 == 0:
                print(f"Epoch {i}, Loss: {loss.item():.4f}")
        
        return self.h_env.detach().cpu().numpy()

# 使用示例:
# solver = scSCF_Solver(cell_line_data, patient_tumor_data, genes)
# solver.define_mask(receptor_gene_indices) # 仅允许受体基因拥有非零场
# fields = solver.train()
# top_fields = fields.argsort()[-10:] # 找到场强最大的前10个受体
