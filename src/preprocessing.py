import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import warnings

# 忽略部分无关紧要的警告
warnings.filterwarnings('ignore')

class scSCF_Preprocessor:
    def __init__(self, vac_path, tissue_path, species='human'):
        """
        vac_path: 细胞系数据路径 (基态, .h5ad)
        tissue_path: 组织/病人数据路径 (受扰态, .h5ad)
        species: 'human' or 'mouse' 用于加载管家基因列表
        """
        print(f"Loading data...\nVacuum (Cell Line): {vac_path}\nTissue (Patient): {tissue_path}")
        self.adata_vac = sc.read_h5ad(vac_path)
        self.adata_tissue = sc.read_h5ad(tissue_path)
        self.species = species
        self.common_genes = []
        
    def qc_and_filter(self, min_genes=200, min_cells=3, mt_cutoff=20):
        """
        基础质控: 过滤低质量细胞和基因
        """
        print("--- Step 1: Quality Control ---")
        for adata, label in zip([self.adata_vac, self.adata_tissue], ['Vacuum', 'Tissue']):
            # 标记线粒体基因
            mt_prefix = 'MT-' if self.species == 'human' else 'Mt-'
            adata.var['mt'] = adata.var_names.str.startswith(mt_prefix)
            sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

            # 过滤
            n_cells_before = adata.n_obs
            sc.pp.filter_cells(adata, min_genes=min_genes)
            sc.pp.filter_genes(adata, min_cells=min_cells)
            adata = adata[adata.obs.pct_counts_mt < mt_cutoff, :]
            
            print(f"[{label}] Cells: {n_cells_before} -> {adata.n_obs}")
            
            # 更新成员变量
            if label == 'Vacuum': self.adata_vac = adata
            else: self.adata_tissue = adata

    def align_features(self):
        """
        取基因交集，确保矩阵维度一致
        """
        print("\n--- Step 2: Feature Alignment (Intersection) ---")
        vac_genes = set(self.adata_vac.var_names)
        tissue_genes = set(self.adata_tissue.var_names)
        self.common_genes = list(vac_genes.intersection(tissue_genes))
        
        # 按照相同顺序截取基因
        self.adata_vac = self.adata_vac[:, self.common_genes].copy()
        self.adata_tissue = self.adata_tissue[:, self.common_genes].copy()
        
        print(f"Aligned Features: {len(self.common_genes)} common genes retained.")

    def robust_normalization(self):
        """
        关键步骤：归一化
        不使用 Harmony/BatchCorrection，而是使用 LogCPM + Housekeeping Scaling
        """
        print("\n--- Step 3: Robust Normalization ---")
        
        # 1. 基础深度归一化 (Target sum = 1e4, 即 CPM/100)
        sc.pp.normalize_total(self.adata_vac, target_sum=1e4)
        sc.pp.normalize_total(self.adata_tissue, target_sum=1e4)
        
        # 2. 对数变换
        sc.pp.log1p(self.adata_vac)
        sc.pp.log1p(self.adata_tissue)
        
        # 3. Housekeeping (HK) Anchor Scaling
        # 加载管家基因列表 (此处模拟一个列表，实际应用建议读取标准 HK 基因集如 Tirosh et al.)
        # 假设我们已经有一个 reliable_hk_genes 列表
        # 在实际操作中，通常使用 ribosomal proteins (RPS/RPL) 或 GAPDH, ACTB 等作为锚点
        hk_genes = [g for g in self.common_genes if g.startswith(('RPS', 'RPL', 'GAPDH', 'ACTB'))]
        
        if len(hk_genes) < 10:
            print("Warning: Too few HK genes found for scaling. Using global mean (risk of bias).")
            hk_genes = self.common_genes # Fallback
            
        print(f"Using {len(hk_genes)} Housekeeping genes for inter-dataset scaling anchors.")
        
        # 计算 HK 基因的平均表达水平
        mean_vac_hk = np.mean(self.adata_vac[:, hk_genes].X, axis=None) # 若是稀疏矩阵需用 .X.mean()
        mean_tissue_hk = np.mean(self.adata_tissue[:, hk_genes].X, axis=None)
        
        # 计算缩放因子 (Scaling Factor)
        # 我们将 Tissue 的强度缩放到与 Vacuum 一致，消除测序平台带来的全局增益差异
        scaling_factor = mean_vac_hk / mean_tissue_hk
        print(f"Scaling Factor (Vacuum/Tissue) = {scaling_factor:.4f}")
        
        # 应用缩放
        # 注意：这里我们线性缩放对数后的值(近似)，或者应该在 log 前缩放 count。
        # 为了物理模型的稳定性，通常在 Z-score 之前做这一步。
        # 这里为了简单，我们直接对 Tissue 数据乘以因子 (假设系统性偏差是乘性的)
        self.adata_tissue.X = self.adata_tissue.X * scaling_factor
        
    def get_matrices_for_solver(self, high_variance_only=True, n_top_genes=2000):
        """
        准备最终输入给 Solver 的 Numpy 矩阵
        """
        print("\n--- Step 4: Exporting Matrices ---")
        
        target_genes = self.common_genes
        
        if high_variance_only:
            # 仅在 Vacuum (基态) 中寻找高变基因，定义系统的“内在自由度”
            sc.pp.highly_variable_genes(self.adata_vac, n_top_genes=n_top_genes, subset=False)
            hvgs = self.adata_vac.var[self.adata_vac.var['highly_variable']].index
            target_genes = list(hvgs)
            print(f"Subsetting to {len(target_genes)} HVGs based on Vacuum state variability.")
        
        # 提取数据
        X_vac = self.adata_vac[:, target_genes].X
        X_tissue = self.adata_tissue[:, target_genes].X
        
        # 处理稀疏矩阵
        if hasattr(X_vac, "toarray"): X_vac = X_vac.toarray()
        if hasattr(X_tissue, "toarray"): X_tissue = X_tissue.toarray()
        
        # 最终中心化 (Z-score centering)
        # 物理模型通常处理相对于平均值的涨落
        # 注意：我们使用 VACUUM 的均值和标准差来标准化 TISSUE
        # 这样 Tissue 的数据就变成了相对于 Vacuum 的“偏移量”
        
        mu_0 = np.mean(X_vac, axis=0)
        std_0 = np.std(X_vac, axis=0) + 1e-6 # 防止除零
        
        X_vac_norm = (X_vac - mu_0) / std_0
        X_tissue_norm = (X_tissue - mu_0) / std_0 # 关键：用基态的参数标准化受扰态
        
        return X_vac_norm, X_tissue_norm, target_genes

# --- 使用示例 ---
# processor = scSCF_Preprocessor('ccle_lung.h5ad', 'patient_lung_tumor.h5ad')
# processor.qc_and_filter()
# processor.align_features()
# processor.robust_normalization()
# X_vac, X_tissue, genes = processor.get_matrices_for_solver()

# 接下来将 X_vac, X_tissue 传入之前的 scSCF_Solver 即可
