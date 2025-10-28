print('Importing packages')

import argparse
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import scanpy as sc
import anndata as ad
import numpy as np
from anndata import AnnData
import pathlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import skimage
import seaborn as sns
import tangram as tg
import logging
from scipy.sparse import issparse
from tangram import utils as ut
import sys
import time

print('Imported packages successfully')

start_time=time.time()

print('Reading in functions plot_auc,construct_obs_plot etc\n')

def random_sampling(adata, n_samples):
    np.random.seed(42)  # For reproducibility
    random_indices = np.random.choice(adata.n_obs, size=n_samples, replace=False)
    return adata[random_indices, :]

def plot_auc(df_all_genes, test_genes=None):
    """
        Plots auc curve which is used to evaluate model performance.
    
    Args:
        df_all_genes (Pandas dataframe): returned by compare_spatial_geneexp(adata_ge, adata_sp); 
        test_genes (list): list of test genes, if not given, test_genes will be set to genes where 'is_training' field is False

    Returns:
        None
    """
    metric_dict, ((pol_xs, pol_ys), (xs, ys)) = ut.eval_metric(df_all_genes, test_genes)
    
    fig = plt.figure()
    plt.figure(figsize=(6, 5))

    plt.plot(pol_xs, pol_ys, c='r')
    sns.scatterplot(xs, ys, alpha=0.5, edgecolors='face')
        
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.gca().set_aspect(.5)
    plt.xlabel('score')
    plt.ylabel('spatial sparsity')
    plt.tick_params(axis='both', labelsize=8)
    plt.title('Prediction on test transcriptome')
    
    textstr = 'auc_score={}'.format(np.round(metric_dict['auc_score'], 3))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.3)
    # place a text box in upper left in axes coords
    plt.text(0.03, 0.1, textstr, fontsize=11, verticalalignment='top', bbox=props);


def construct_obs_plot(df_plot, adata, perc=0, suffix=None):
    # clip
    df_plot = df_plot.clip(df_plot.quantile(perc), df_plot.quantile(1 - perc), axis=1)

    # normalize
    df_plot = (df_plot - df_plot.min()) / (df_plot.max() - df_plot.min())

    if suffix:
        df_plot = df_plot.add_suffix(" ({})".format(suffix))
    adata.obs = pd.concat([adata.obs, df_plot], axis=1)

def plot_cell_annotation_sc(
    adata_sp, 
    annotation_list, 
    x="center_x", 
    y="center_y", 
    spot_size=None, 
    scale_factor=None, 
    perc=0,
    alpha_img=1.0,
    bw=False,
    ax=None
):

    # Ensure only the first annotation is used
    if isinstance(annotation_list, list):
        annotation = annotation_list[0]
    else:
        annotation = annotation_list
        
    # remove previous df_plot in obs
    adata_sp.obs.drop([annotation], inplace=True, errors="ignore", axis=1)

    # construct df_plot
    df = adata_sp.obsm["tangram_ct_pred"][[annotation]]
    construct_obs_plot(df, adata_sp, perc=perc)
    
    #non visium data 
    if 'spatial' not in adata_sp.obsm.keys():
        #add spatial coordinates to obsm of spatial data 
        coords = [[x,y] for x,y in zip(adata_sp.obs[x].values,adata_sp.obs[y].values)]
        adata_sp.obsm['spatial'] = np.array(coords)
        
    # Check for spot size and scale factor if spatial data is not present
    if 'spatial' not in adata_sp.uns.keys():
        if spot_size is None or scale_factor is None:
            raise ValueError("Spot Size and Scale Factor cannot be None when adata_sp.uns['spatial'] does not exist")
    else:
        # Ensure spot size and scale factor are None if spatial data is present
        if spot_size is not None or scale_factor is not None:
            raise ValueError("Spot Size and Scale Factor should be None when adata_sp.uns['spatial'] exists")

    sc.pl.spatial(
        adata_sp, color=[annotation], cmap="viridis", show=False, frameon=False, spot_size=spot_size,
        scale_factor=scale_factor, alpha_img=alpha_img, bw=bw, ax=ax
    )

    adata_sp.obs.drop([annotation], inplace=True, errors="ignore", axis=1)

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Get cell type annotation for spatial transcriptomics data')
parser.add_argument('--path', required=True, help='Path to the sc and st data')
parser.add_argument('--prefix', required=True, help='Prefix to the cell and gene metadata for spatial')
parser.add_argument('--file_sc', required=True, help='Filename of the single cell data')
parser.add_argument('--n_samples_sc', type=int, default=30000, help='Number of single cells to consider')
parser.add_argument('--start_cell', type=int, default=0, help='Starting cell index for single cell data')
parser.add_argument('--end_cell', type=int, help='Ending cell index for single cell data')

#parser.add_argument('--n_samples_st', type=int, default=5000, help='Number of spatial cells to consider')

#n_samples_sc

args = parser.parse_args()

path = args.path
prefix = args.prefix
file_sc = args.file_sc
#file_st = args.file_st
n_samples_sc = args.n_samples_sc
#n_samples_st = args.n_samples_st
start_cell=args.start_cell
end_cell=args.end_cell

print('Loading single cell cancer data' + '/path/to/file/'+file_sc)

adata_sc = ad.read_h5ad('/path/to/file/'+file_sc)

print('Read in adata_sc successfully')

print('adata_sc dimensions:', adata_sc.shape)
if adata_sc.shape[0] == 0 or adata_sc.shape[1] == 0:
    raise ValueError("adata_sc is empty. Please check the input data.")

print('adata_sc:\n')
print(adata_sc)

print('adata_sc.obs[annotation_major]:\n')
print(adata_sc.obs['annotation_major'])

print('np.unique(adata_sc.obs[annotation_major])')
print(np.unique(adata_sc.obs['annotation_major']))

cell_by_gene = pd.read_csv(path+prefix+'cell_by_gene.csv')
cell_metadata = pd.read_csv(path+prefix+'cell_metadata.csv')

adata_st = ad.AnnData(cell_by_gene)

print('Getting raw st data')

raw_data=adata_st.X.copy()

adata_st.layers['counts']=raw_data

print('Here is the spatial raw count data:\n')
print(adata_st.layers['counts'])

adata_st.obs=cell_metadata

print('adata_sc has {} cells and {} genes'.format(adata_sc.n_obs, adata_sc.n_vars))
print('adata_st has {} cells and {} genes'.format(adata_st.n_obs, adata_st.n_vars))

adata_sc = adata_sc[:n_samples_sc, :]

print('start_cell')
print(start_cell)
print('end_cell')
print(end_cell)

# Slice the single cell data based on start_cell and end_cell
if end_cell is not None:
    adata_st = adata_st[start_cell:end_cell, :]
#else:
#    adata_st = adata_st[start_cell:start_cell + n_samples_sc, :]

#adata_st = adata_st[:n_samples_st, :]

print('adata_sc has {} cells and {} genes after subsetting according to n_samples_sc'.format(adata_sc.n_obs, adata_sc.n_vars))
print('adata_st has {} cells and {} genes after subsetting according to n_samples_st'.format(adata_st.n_obs, adata_st.n_vars))

print('Ensure raw attribute of sc is not causing issues')
if adata_sc.raw is not None:
    print("Clearing raw attribute from adata_sc")
    adata_sc.raw = None

#print('Filtering cells, min_genes=20 for both sc and st')
print('Running sc.pp.filter_cells(adata_st, min_genes=10)')
sc.pp.filter_cells(adata_st, min_genes=10)

print('Running sc.pp.filter_cells(adata_sc, min_genes=10)')
sc.pp.filter_cells(adata_sc, min_genes=10)

print('Filtered adata_sc has {} cells and {} genes'.format(adata_sc.n_obs, adata_sc.n_vars))

print('Filtered adata_st has {} cells and {} genes'.format(adata_st.n_obs, adata_st.n_vars))

print('Filtering genes, min_cells=3 for both sc and st')

sc.pp.filter_genes(adata_sc, min_cells=3)
sc.pp.filter_genes(adata_st, min_cells=3)

print('Filtered adata_sc has {} cells and {} genes'.format(adata_sc.n_obs, adata_sc.n_vars))

print('Filtered adata_st has {} cells and {} genes'.format(adata_st.n_obs, adata_st.n_vars))

print('Normalising and log1p for adata_st: \n')

sc.pp.normalize_total(adata_st)
sc.pp.log1p(adata_st)

print('Here is the normalised st data:\n')

print(adata_st.X)

print('Looking at adata_sc.X:\n')
print(adata_sc.X)

print('Data already appears to be log normalised')

print('Selecting the top 2000 highly variable genes in sc')
sc.pp.highly_variable_genes(adata_sc, n_top_genes=2000)
adata_sc = adata_sc[:, adata_sc.var.highly_variable]

print('Filtered adata_sc has {} cells and {} genes'.format(adata_sc.n_obs, adata_sc.n_vars))
print('adata_st has {} cells and {} genes'.format(adata_st.n_obs, adata_st.n_vars))

print('Converting to float32 data type')
adata_sc.X = adata_sc.X.astype(np.float32)
adata_st.X = adata_st.X.astype(np.float32)

# Plotting
print('Generating and saving plots')
fig, ax = plt.subplots(1, 1, figsize=(10, 5))

# Plot UMAP
sc.pl.umap(
    adata_sc, color="annotation_major", size=10, frameon=False, show=False, ax=ax
)

plt.tight_layout()
plot_file = path + 'umap_sc_plot_annotation_major.png'
plt.savefig(plot_file)
print(f'Plots saved to {plot_file}')

print('Running rank_genes_groups')

sc.tl.rank_genes_groups(adata_sc, groupby="annotation_major", use_raw=False)

print('Getting markers')

#exit()

markers_df = pd.DataFrame(adata_sc.uns["rank_genes_groups"]["names"]).iloc[0:100, :]
markers = list(np.unique(markers_df.melt().value.values))

print('Length of markers:\n')
print(len(markers))

#print('Exiting!')
#exit()

print('Running tg.pp_adatas(adata_sc, adata, genes=markers)')

tg.pp_adatas(adata_sc, adata_st, genes=markers)

print('Running tg.map_cells_to_space with 600 epochs')

ad_map = tg.map_cells_to_space(adata_sc, adata_st,
    mode="cells",
    density_prior='uniform',
    num_epochs=600,
    device='cpu',
)

print('Done with tg.map_cells_to_space!\n')

print('Here is ad_map:\n')
print(ad_map)

print('Here is adata_st before projecting the cell annotations:\n')
print(adata_st)

print('\n')
print('Project cell annotations to space')
tg.project_cell_annotations(ad_map, adata_st, annotation="annotation_major")
print('Projected cell annotations to spatial data successfully')

print('\n')

print('Here is adata_st after projecting the cell annotations:\n')
print(adata_st)

# Verify the projection
print('Updated adata_st.obs with projected cell annotations:\n')
#print(adata_st.obs[['annotation_major']].head())

print(adata_st.obsm['tangram_ct_pred'])

#print(np.unique(adata_st.obsm['tangram_ct_pred']))

# Get the probabilities matrix
tangram_ct_pred = adata_st.obsm['tangram_ct_pred']

# Find the index of the maximum probability for each row
max_prob_indices = np.argmax(tangram_ct_pred, axis=1)

# Map the indices to cell type names
cell_types = tangram_ct_pred.columns[max_prob_indices]

# Add the predicted cell types to the obs dataframe
adata_st.obs['tangram_ct'] = cell_types

print('Writing split h5ad file')
adata_st.write_h5ad(path+prefix+'adata_st_'+str(start_cell)+'_split.h5ad')

print('Exiting!')
exit()
