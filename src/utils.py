import numpy as np
import torch
import pandas as pd
import matplotlib.pyplot as plt

def generate_robust_synthetic_data(n_genes=200, n_vac=800, n_tissue=3000, n_active=5):
    """
    生成仿真数据：返回活跃受体索引，模拟不对称样本量和响应异质性
    """
    # 1. 构造基态协方差 C0 (代表基因内在耦合)
    A = np.random.randn(n_genes, n_genes) * 0.1
    C0 = np.dot(A, A.T) + np.eye(n_genes) * 0.5
    
    # 2. 生成基态细胞 (Vacuum)
    X_vac = np.random.multivariate_normal(np.zeros(n_genes), C0, size=n_vac)
    
    # 3. 设定 Ground Truth 场
    h_gt = np.zeros(n_genes)
    active_indices = np.random.choice(n_genes, n_active, replace=False)
    # 给激活的受体赋予随机正负场强（模拟激活或抑制）
    h_gt[active_indices] = np.random.uniform(low=2.0, high=5.0, size=n_active) * np.random.choice([-1, 1], n_active)
    
    # 4. 计算受扰态的理论偏移: delta_mu = C0 @ h
    delta_mu = np.dot(C0, h_gt)
    
    # 5. 生成组织细胞 (Tissue) - 模拟 70% 的响应率
    n_responders = int(n_tissue * 0.7)
    n_non_responders = n_tissue - n_responders
    
    # 响应群体（均值偏移）
    X_resp = np.random.multivariate_normal(delta_mu, C0, size=n_responders)
    # 非响应群体（均值不偏移，模拟远离互作现场的细胞）
    X_non_resp = np.random.multivariate_normal(np.zeros(n_genes), C0, size=n_non_responders)
    
    X_tissue = np.vstack([X_resp, X_non_resp])
    
    print(f"--- Simulation Config ---")
    print(f"Vacuum Cells: {n_vac}, Tissue Cells: {n_tissue}")
    print(f"Response Rate: 70%, Active Receptors: {n_active}")
    
    return X_vac, X_tissue, h_gt, active_indices
