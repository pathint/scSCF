# 生成测试数据
X_vac, X_tissue, h_gt, active_idx = generate_robust_synthetic_data()

# 初始化 Solver 并训练
# 假设我们知道潜在受体列表（实际研究中从数据库获取，这里随机选50个包含正确答案的基因）
candidate_receptors = list(set(list(active_idx) + list(np.random.choice(200, 45, replace=False))))
solver = scSCF_Solver(X_vac, X_tissue, [f"Gene_{i}" for i in range(200)])
solver.define_mask(candidate_receptors)
h_pred = solver.train(epochs=600, lr=0.05, l1_lambda=0.3)

# --- 验证结果可视化 ---
plt.figure(figsize=(10, 4))

# 子图1: 场强的回归对比
plt.subplot(1, 2, 1)
plt.scatter(h_gt, h_pred, alpha=0.6, color='teal')
plt.xlabel("Ground Truth Field ($h_{gt}$)")
plt.ylabel("Predicted Field ($h_{pred}$)")
plt.title("Field Recovery Accuracy")
plt.plot([0, 5], [0, 5], '--', color='red') # 对角线

# 子图2: Top 基因识别情况
plt.subplot(1, 2, 2)
df_res = pd.DataFrame({'GT': h_gt, 'Pred': h_pred})
df_res = df_res.sort_values('Pred', ascending=False).head(20)
colors = ['red' if x > 0 else 'gray' for x in df_res['GT']]
plt.bar(range(20), df_res['Pred'], color=colors)
plt.title("Top 20 Predicted Interaction Channels\n(Red = True Positives)")
plt.xlabel("Rank")
plt.ylabel("Field Intensity")

plt.tight_layout()
plt.show()
