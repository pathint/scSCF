# scSCF：single-cell RNA self-consistent field interaction solver 

Self-consistent field solver for cell - signaling field interactions based on single-cell RNA sequencing profiles

## 算法框架 （Algorithm Flowchart）

#### 1. 核心物理假设与映射（Assumputions）

我们将肿瘤细胞在微环境中的基因表达视为一个受“外场”影响的系统。

* **粒子（Particle）：** 单个肿瘤细胞。
* **内能（Internal Hamiltonian, $H_0$）：** 肿瘤细胞内在的基因调控网络（GRN）。这可以通过纯肿瘤细胞系（Cell Line）的数据来表征（即不受微环境干扰的“真空”状态）。
* **外场（External Field, $\Phi$）：** 由肿瘤微环境（TME）中其他细胞（免疫细胞、基质细胞）分泌的配体（Ligand）形成的“平均场”。
* **相互作用能（Interaction Energy, $H_{int}$）：** 肿瘤细胞受体（Receptor）与外场配体结合后，引发下游信号通路改变表达谱的能量消耗。

#### 2. 数据输入 (Inputs)

1. **$D_{vac}$ (Vacuum State):** 肿瘤细胞系单细胞数据（无微环境，仅含内在调控）。
2. **$D_{env}$ (Environmental State):** 组织或类器官来源的单细胞数据（包含肿瘤细胞及微环境细胞）。
3. **$DB_{LR}$:** 配体-受体相互作用数据库（如CellChatDB, CellPhoneDB）。

#### 3. 算法流程 (Workflow)

算法的核心是通过SCF思路，重构肿瘤细胞在组织中的表达谱。如果我们可以用“细胞系表达 + TME场效应”完美重构出“组织内肿瘤细胞表达”，那么这个“TME场效应”就是我们要找的关键因素。

#### 3.1：构建基态哈密顿量 ($H_0$)

利用 $D_{vac}$（细胞系数据），构建肿瘤细胞的内在基因共表达概率分布。基于最大熵原理（MaxEnt），细胞处于某一基因表达状态 $\mathbf{x}$ 的概率为：

$$P_0(\mathbf{x}) = \frac{1}{Z_0} \exp(-H_0(\mathbf{x}))$$

其中， $H_0$ 可采用类似 Ising 模型或 Gaussian Graphical Model 的形式，描述基因间的内在耦合：

$$H_0(\mathbf{x}) = -\sum_i h_i x_i - \sum J_{ij} x_i x_j$$

* $x_i$: 基因 $i$ 的表达量。
* $J_{ij}$: 基因 $i$ 和 $j$ 的内在共表达强度（从细胞系数据训练得到）。

#### 3.2：定义平均场与微环境微扰 ($H_{eff}$)

在组织数据 $D_{env}$ 中，肿瘤细胞受到微环境细胞 $k$ (如T细胞, 成纤维细胞) 的影响。将这种影响定义为外场微扰。

有效哈密顿量变为：

$$H_{eff}(\mathbf{x}) = H_0(\mathbf{x}) - \sum_{r \in Receptors} \lambda_r \cdot \phi_r \cdot x_r$$

* $x_r$: 肿瘤细胞中受体基因 $r$ 的表达量（或其下游靶基因的集合）。
* $\phi_r$: 平均场强度。这是我们试图求解的核心。它代表了环境中有多少配体能够激活受体 $r$。
* $\phi_r \approx \sum_{k \in CellTypes} w_k \cdot L_{k,r}$
* $L_{k,r}$: 细胞类型 $k$ 表达配体的平均水平。
* $w_k$: 细胞类型 $k$ 在微环境中的密度/空间邻近度。
* $\lambda_r$: 耦合常数（Susceptibility），表示受体 $r$ 被激活后对整个转录组 $\mathbf{x}$ 的扰动能力。

#### 3.3：自洽场迭代 (SCF Iteration)

我们的目标是找到一组最优的场参数 $\{\phi_r\}$ 和响应系数 $\{\lambda_r\}$，使得重构出的分布 $P_{eff}$ 最接近真实的组织内肿瘤细胞分布 $P_{tissue}$。

1. **初始化：** 设定初始场 $\Phi^{(0)}$ 为 0 或随机小量。
2. **预测 (Prediction):** 基于当前场 $\Phi^{(t)}$，计算理论上的肿瘤细胞表达谱期望值：

$$\langle \mathbf{x} \rangle_{model} = \sum \mathbf{x} \cdot P(\mathbf{x} | H_0, \Phi^{(t)})$$

这一步可以通过MCMC采样或变分推断近似计算。

3. **校正 (Correction/Self-Consistency):** 比较预测表达谱 $\langle \mathbf{x} \rangle_{model}$ 与真实的组织来源肿瘤细胞表达谱 $\langle \mathbf{x} \rangle_{tissue}$。
   
* 计算误差函数（如 KL 散度或欧氏距离）：

$$\mathcal{L} = || \langle \mathbf{x} \rangle_{model} - \langle \mathbf{x} \rangle_{tissue} ||^2$$

4. **更新场 (Update):** 利用梯度下降法更新场强度，使其满足自洽条件：

$$\Phi^{(t+1)} = \Phi^{(t)} - \eta \frac{\partial \mathcal{L}}{\partial \Phi}$$

5. **收敛：** 当 $\Delta \Phi$ 小于阈值时停止。此时的 $\Phi^*$ 即为**重构出的TME关键互作场**。

#### 3.4：关键因素解析 (Decoupling)

算法收敛后，我们得到最终的 $\Phi^*_r$（针对每个受体的有效场强），可以：

* **识别关键受体：** $\Phi^*_r$ 值最大的受体即为介导TME影响的关键通道。
* **溯源关键细胞：** 将 $\Phi^*_r$ 分解回细胞来源 $\sum w_k L_{k,r}$。
* 如果 $\Phi^*_r$ 很高，且主要由成纤维细胞的配体贡献，则推断**成纤维细胞-肿瘤细胞互作**是导致肿瘤从“细胞系状态”向“组织状态”转变的主因。

