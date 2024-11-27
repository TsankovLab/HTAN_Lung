import os, glob, re, pickle
from functools import partial
from collections import OrderedDict
import operator as op
from cytoolz import compose

import pandas as pd
import seaborn as sns
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt

from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_binarization, plot_rss

from IPython.display import HTML, display

DATA_INPUT = "./data/"
FIGURES_FOLDERNAME = "../../figures/Fig3_stroma/"

sc.settings.figdir = FIGURES_FOLDERNAME

#####
auc_mtx = pd.read_csv(os.path.join(DATA_INPUT, "HTAN.auc.csv"), index_col=0)

mes = sc.read(os.path.join(DATA_INPUT, "mes.subset.bysubtype.h5ad"))
mes.obs

rss = regulon_specificity_scores(auc_mtx, mes.obs.subtype)
rss

##### Extended 7D
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 3), dpi=100)
plot_rss(rss, 'CAF.COL', ax=ax1)
ax1.set_xlabel('')
plot_rss(rss, 'CAF.ISG', ax=ax2)
ax2.set_ylabel('')
plot_rss(rss, 'CAF.Cycling', ax=ax3)
ax3.set_xlabel('')
ax3.set_ylabel('')

plt.savefig(os.path.join(FIGURES_FOLDERNAME, "Ext_Fig7D.pdf"), format='pdf', dpi=100)