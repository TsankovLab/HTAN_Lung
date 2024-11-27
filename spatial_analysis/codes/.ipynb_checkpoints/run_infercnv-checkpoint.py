import os
import anndata as ad
import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt

from tqdm import tqdm
from pathlib import Path
from scanpy._settings import ScanpyConfig

plt.rcParams['figure.dpi'] = 200

# init constants
FEATURE = 'celltype'

slides = ['BW09_D', 'BW11_C', 'BW14_B', 'MGH1174_C', 'MGH1174_D']

OUT_DIR = "./data/infercnv_input"
FIG_DIR = "../../figures/cnv_figures"
ScanpyConfig.figdir = Path(FIG_DIR)

merged_adata = None

for slide in tqdm(slides):
    print()
    print('#'*15)
    print(f'Processing {slide}...')
    print('#'*15)
    print()

    adata_st = ad.read_h5ad(os.path.join(OUT_DIR, f'{slide}.h5ad'))

    sc.pp.calculate_qc_metrics(adata_st, inplace=True)
    
    if slide == 'BW09_D':
        # filter bottom rows - bad quality (only for BW09.D)
        adata_st = adata_st[(adata_st.obs['array_row']>15)].copy()

    adata_st.obs['sample_id'] = slide
    
    cnv.io.genomic_position_from_gtf(os.path.join(OUT_DIR, 'ref', 'Homo_sapiens.GRCh38.111.gtf'), adata_st, gtf_gene_id='gene_id', adata_gene_id='gene_ids')

    no_rm = [gene for gene in adata_st.var.index if gene not in adata_st.var.loc[adata_st.var['start'].isna()].index.tolist()]
    adata_st = adata_st[:, no_rm].copy()

    adata_st.var['start'] = adata_st.var['start'].astype(int)
    adata_st.var['end'] = adata_st.var['end'].astype(int)

    non_tumor_ctypes = adata_st.obs[f'{FEATURE}_predicted'].unique().tolist()
    non_tumor_ctypes.remove('Cancer')

    adata_st.var['chromosome'] = list(map(lambda x: f'chr{x}', adata_st.var['chromosome']))

    cnv.tl.infercnv(adata_st, window_size=250, layer='raw', exclude_chromosomes=adata_st.var['chromosome'].unique().tolist()[22:])

    cnv.tl.pca(adata_st)
    cnv.pp.neighbors(adata_st)
    cnv.tl.leiden(adata_st)
    cnv.tl.umap(adata_st)
    cnv.tl.cnv_score(adata_st)

    sc.pl.spatial(adata_st, color='cnv_score', alpha_img=0.6, vmax=0.02, show=False, save=f'{slide}_cnv_score.pdf')

print('DONE!')
