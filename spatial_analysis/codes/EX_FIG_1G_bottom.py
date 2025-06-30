import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

plt.rcParams['figure.dpi'] = 350

OUT_DIR = "./data/infercnv_input"
FIG_DIR = "../../figures/"
os.makedirs(FIG_DIR, exist_ok=True)

merged_df = pd.read_csv('EX_FIG_1E_bottom.csv')

dataset_order = ['Zhao et al. this study', 'Zuani et al. Nat Comm 2024', 'Takano et al. Nat Comm 2024', 'Li et al. Cancer Cell 2022', 'Greenwald wt al. Cell 2024']
cmap = {
    'Zhao et al. this study': 'black',
    'Zuani et al. Nat Comm 2024': 'tab:gray',
    'Takano et al. Nat Comm 2024': 'tab:gray',
    'Li et al. Cancer Cell 2022': 'tab:gray',
    'Greenwald wt al. Cell 2024': 'tab:gray',
}

fig, ax = plt.subplots(1, 1, figsize=(1,1))
sns.boxplot(merged_df, x='dataset', y='log10_nCount_Spatial', hue='dataset', palette=cmap, order=dataset_order, showcaps=False, flierprops={"marker": ""}, linewidth=1.1, fill=False, width=0.5)
plt.xlabel('Dataset', fontsize=6)
plt.ylabel('Number of features\nper spot', fontsize=6)
plt.xticks(ticks=[0, 1, 2, 3, 4], labels=dataset_order, rotation=45, ha='right', fontsize=6)
plt.yticks(ticks=[0, 1, 2, 3, 4, 5], labels=[0, 10, 100, 1000, 10000, 100000], fontsize=6)
ax.legend().remove()
plt.savefig(os.path.join(FIG_DIR, 'EX_FIG_1G_bottom.pdf'), bbox_inches='tight')
plt.close()