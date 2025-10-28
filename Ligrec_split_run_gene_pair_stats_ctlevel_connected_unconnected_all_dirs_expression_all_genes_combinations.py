print('Importing packages')

import os
import numpy as np
import pandas as pd
import itertools
from itertools import permutations
import scanpy as sc
#import squidpy as sq
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import anndata as ad
from pathlib import Path
import scipy
from matplotlib.pyplot import rc_context
import os
import csv
sc.logging.print_header()
import random
import argparse
from scipy.stats import pearsonr, spearmanr
from scipy import stats
import tracemalloc
import time
#sc.logging.print_header()

# Record start time
start_time = time.time()

# Create the argument parser
parser = argparse.ArgumentParser(description='Find genes expressed in cancer cells connected to normal cells, and vice versa, when splitting up a large list of gene combinations')

# Add the --path argument
parser.add_argument('--path', type=str, help='Path to files')
parser.add_argument('--cell_subtype', type=str, help='Normal cell subtype')
parser.add_argument('--split_part', type=int, required=True, help='Split part (1, 2, 3, ...)')
parser.add_argument('--split_size', type=int, required=True, help='Size of split e.g. how many gene-pairs to consider in a chunk')
parser.add_argument('--start_index', type=int, required=True, help='Which index to start on across the gene combinations')
parser.add_argument('--end_index', type=int, required=True, help='Which index to end on across the gene combinations')
parser.add_argument('--ct_level', type=str, required=True, help='Cell type annotation')
parser.add_argument('--version', type=str, required=True, help='Run version')
parser.add_argument('--corr_coef', type=str, choices=['pearson', 'spearman'], default='spearman', help='Correlation coefficient method (pearson or spearman)')
parser.add_argument('--sample',type=str,help="Which sample")
parser.add_argument('--all_genes', action='store_true', help='Consider all genes')
parser.add_argument('--ligrec_known', action='store_true', help='Consider known lig-rec pairs')
parser.add_argument('--DE_known', action='store_true', help='Consider DE gene-pairs')

args = parser.parse_args()

if args.all_genes:
    print("Considering all gene pairs")
    all_genes=True
if args.ligrec_known:
    print("Considering known ligrec_pairs")
    ligrec_known=True
if args.DE_known:
    print("Considering DE-genes")
    DE_known=True

# Access the value of the --path argument
path = args.path
cell_subtype = args.cell_subtype
split_part=args.split_part
split_size=args.split_size
start_index=args.start_index
end_index=args.end_index
ct_level=args.ct_level
version=args.version
sample=args.sample
all_genes=args.all_genes
ligrec_known=args.ligrec_known
DE_known=args.DE_known
corr_coef=args.corr_coef

print('Printing command line arguments')

print('path:\n')
print(path)
print('cell_subtype:\n')
print(cell_subtype)
print('split part:\n')
print(split_part)
print('split size:\n')
print(split_size)
print('start_index:\n')
print(start_index)
print('end index:\n')
print(end_index)
print('ct_level:\n')
print(ct_level)
print('version:\n')
print(version)
print('all_genes:\n')
print(all_genes)
print('ligrec_known:\n')
print(ligrec_known)
print('DE_known:\n')
print(DE_known)
print('corr_coef:\n')
print(corr_coef)
print('sample:\n')
print(sample)

if (all_genes):
    all_flag='all_genes'
if (ligrec_known):
    all_flag='known_ligrec_genes'
if (DE_known):
    all_flag='DE_genes'

print('all_flag')
print(all_flag)

if corr_coef == 'pearson':
    # Use Pearson correlation coefficient
    correlation_function = pearsonr
elif corr_coef == 'spearman':
    # Use Spearman correlation coefficient
    correlation_function = spearmanr

def lig_rec_plot_normal_cancer_CN(adata,cell_subtype,impute=True,lig_rec_pair=('PDGFB', 'LRP1'),all_genes=False,corr_coef=spearmanr):
    if (impute):
        #print('Using imputed values')
        adata.X=adata.obsm['X_normalized_scVI']
        assert(np.sum(np.asarray(adata.X))==np.sum(np.asarray(adata.obsm['X_normalized_scVI'])))
    else:
        print('Not using imputed values')
    print('Considering '+Normal_cells[cell_subtype])
    print('Cancer_cell')
    print(Cancer_cell)
    adata_cancer_subtype=adata[(adata.obs[ct_level]==(Cancer_cell)) | (adata.obs[ct_level]==(Normal_cells[cell_subtype]))]
    print('adata_cancer_subtype.obsp[spatial_connectivities]:\n')
    print(adata_cancer_subtype.obsp['spatial_connectivities'])
    # Check if 'connectivities' exists in obsp and if it's empty
    if 'connectivities' in adata_cancer_subtype.obsp:
        connectivities_empty = adata_cancer_subtype.obsp['connectivities'].nnz == 0  # Check if the number of non-zero elements is 0
    else:
        connectivities_empty = True
    # If 'connectivities' is empty, set it to 'spatial_connectivities'
    if connectivities_empty:
        print('Note that connectivities was empty, setting it to spatial_connectivities')
        adata_cancer_subtype.obsp['connectivities'] = adata_cancer_subtype.obsp['spatial_connectivities']
    assert(len(adata_cancer_subtype)==adata_cancer_subtype.obsp['connectivities'].shape[0])
    assert(len(adata_cancer_subtype)==adata_cancer_subtype.obsp['connectivities'].shape[1])
    adjacency_matrix=adata_cancer_subtype.obsp['connectivities'].toarray()
    adjacency_matrix_df=pd.DataFrame(adjacency_matrix,index=adata_cancer_subtype.obs.index, columns=adata_cancer_subtype.obs.index)
    index_list1 = adata_cancer_subtype.obs[(adata_cancer_subtype.obs[ct_level] == Cancer_cell)].index
    index_list2 = adata_cancer_subtype.obs[(adata_cancer_subtype.obs[ct_level]==Normal_cells[cell_subtype])].index
    index_list1 = index_list1.astype(str)
    index_list2 = index_list2.astype(str)
    adjacency_matrix_rows = adjacency_matrix_df.loc[index_list2, index_list1]
    row_labels = adjacency_matrix_rows.index[np.where(adjacency_matrix_rows == 1.0)[0]]
    col_labels = adjacency_matrix_rows.columns[np.where(adjacency_matrix_rows == 1.0)[1]]
    row_labels_uc = adjacency_matrix_rows.index[np.where((adjacency_matrix_rows != 1.0).all(axis=1) & (adjacency_matrix_rows.sum(axis=1) == 0))[0]]
    col_labels_uc = adjacency_matrix_rows.columns[np.where((adjacency_matrix_rows != 1.0).all(axis=0) & (adjacency_matrix_rows.sum(axis=0) == 0))[0]]
    ### For CN direction
    cancer_ligand=lig_rec_pair[0]
    normal_receptor=lig_rec_pair[1]
    idx_ligand=adata.var_names.str.upper().get_loc(cancer_ligand)
    idx_receptor=adata.var_names.str.upper().get_loc(normal_receptor)
    cancer_ligand_list = adata_cancer_subtype[col_labels, :].X[:, idx_ligand]
    normal_receptor_list = adata_cancer_subtype[row_labels, :].X[:, idx_receptor]
    mean_cancer_list=np.mean(cancer_ligand_list)
    mean_normal_list=np.mean(normal_receptor_list)
    corr_connected=round(corr_coef(cancer_ligand_list, normal_receptor_list)[0], 4)
    #print('corr_connected CN:',corr_connected)
    sig_corr_connected=round(corr_coef(cancer_ligand_list, normal_receptor_list)[1], 4)
    connected_df_CN=pd.DataFrame()
    connected_df_CN['normal_receptor_list']=normal_receptor_list
    connected_df_CN['cancer_ligand_list']=cancer_ligand_list
    n_con_cancer=len(cancer_ligand_list)
    n_con_normal=len(normal_receptor_list)
    mean_cancer_ligand=np.mean(cancer_ligand_list)
    mean_normal_receptor=np.mean(normal_receptor_list)
    max_cancer_ligand=np.max(cancer_ligand_list)
    max_normal_receptor=np.max(normal_receptor_list)
    min_cancer_ligand=np.min(cancer_ligand_list)
    min_normal_receptor=np.min(normal_receptor_list)
    cancer_ligand_list_uc = adata_cancer_subtype[col_labels_uc,:].X[:,idx_ligand]
    normal_receptor_list_uc = adata_cancer_subtype[row_labels_uc,:].X[:,idx_receptor]
    cancer_ligand_list_uc_repeats=np.repeat(cancer_ligand_list_uc,len(normal_receptor_list_uc))
    normal_receptor_list_uc_repeats=np.repeat(normal_receptor_list_uc,len(cancer_ligand_list_uc))
    normal_receptor_list_uc = np.array(normal_receptor_list_uc)
    cancer_ligand_list_uc = np.array(cancer_ligand_list_uc)
    N = 0.1*len(normal_receptor_list_uc)*len(cancer_ligand_list_uc)
    sampled_indices1 = np.random.choice(len(cancer_ligand_list_uc), int(N))
    sampled_indices2 = np.random.choice(len(normal_receptor_list_uc), int(N))
    sampled_cancer_ligand_list_uc = cancer_ligand_list_uc[sampled_indices1]
    sampled_normal_receptor_list_uc = normal_receptor_list_uc[sampled_indices2]
    unconnected_df_CN=pd.DataFrame()
    unconnected_df_CN['normal_receptor_list_uc']=sampled_normal_receptor_list_uc
    unconnected_df_CN['cancer_ligand_list_uc']=sampled_cancer_ligand_list_uc
    n_uncon_cancer=len(sampled_cancer_ligand_list_uc)
    n_uncon_normal=len(sampled_normal_receptor_list_uc)
    corr_unconnected=round(corr_coef(sampled_cancer_ligand_list_uc, sampled_normal_receptor_list_uc)[0], 4)
    sig_corr_unconnected=round(corr_coef(sampled_cancer_ligand_list_uc, sampled_normal_receptor_list_uc)[1], 4)
    if ((corr_connected > 0.05) or (corr_connected < -0.05)):
        print('Plotting scatterplot for potential ligand-receptor pair')
        plt.close()
        plt.figure(figsize=(10, 5))
        plt.subplot(1, 2, 1)
        plt.violinplot(cancer_ligand_list,showmedians=True)
            #plt.text(0.5, -0.4, lig_rec_pair[0], horizontalalignment='center', verticalalignment='center', fontsize=12)
        plt.text(0.5, min(cancer_ligand_list) - 0.1, lig_rec_pair[0], horizontalalignment='center', verticalalignment='center', fontsize=12)
        plt.subplot(1, 2, 2)
        plt.violinplot(normal_receptor_list,showmedians=True)
        plt.text(0.5, min(normal_receptor_list) - 0.1, lig_rec_pair[1], horizontalalignment='center', verticalalignment='center', fontsize=12)
        plt.savefig(path+Normal_cells[cell_subtype]+'_imgs/'+'Violin_plot_'+Normal_cells[cell_subtype]+'_'+ct_level+ '_'+ version+'_'+lig_rec_pair[0]+'-'+lig_rec_pair[1]+'.png')
        plt.close()
        plt.hist2d(normal_receptor_list, cancer_ligand_list, bins=30, cmap='viridis')
        plt.colorbar(label='Frequency')
        plt.xlabel(lig_rec_pair[1])
        plt.ylabel(lig_rec_pair[0])
        plt.title('2D Histogram of ligand-receptor expression')
        plt.savefig(path+Normal_cells[cell_subtype]+'_imgs/'+'2D_hist_'+Normal_cells[cell_subtype]+'_'+ct_level+ '_'+ version+'_'+lig_rec_pair[0]+'-'+lig_rec_pair[1]+'.png')
        plt.close()
        plt.figure(figsize=(10, 5))
        try:
            sns.kdeplot(data=connected_df_CN, color="blue")
        except ValueError as e:
            print("Error occurred while plotting KDE:", e)
            print("Gene pair responsible is:\n")
            print(lig_rec_pair[0]+' - '+lig_rec_pair[1])
        #plt.close()
        plt.title('C->'+Normal_cells[cell_subtype]+' '+lig_rec_pair[0]+' - '+lig_rec_pair[1])
        plt.scatter(normal_receptor_list,cancer_ligand_list,color='blue',label='Connected pairs, correlation= '+str(corr_connected), alpha=0.5)
        plt.xlabel(Normal_cells[cell_subtype]+ ' mean = '+str(np.mean(normal_receptor_list)))
        plt.ylabel('Epithelial '+' mean = '+str(np.mean(cancer_ligand_list)))
        plt.legend()
        if (all_genes):
            print('Doing all_genes!') 
            if ((corr_connected>0.8) or (corr_connected < -0.8)):
                plt.savefig(path+Normal_cells[cell_subtype]+'_imgs/'+'CN_Connected_pair_'+Normal_cells[cell_subtype]+'_'+ct_level+ '_'+ version+'_'+lig_rec_pair[0]+'-'+lig_rec_pair[1]+'_cancer_ligand.png')
        else:
            plt.savefig(path+Normal_cells[cell_subtype]+'_imgs/'+'CN_Connected_pair_'+Normal_cells[cell_subtype]+'_'+ct_level+ '_'+ version+'_'+lig_rec_pair[0]+'-'+lig_rec_pair[1]+'_cancer_ligand.png')
    return {'corr_connected':corr_connected, 'sig_corr_connected':sig_corr_connected, 'corr_unconnected':corr_unconnected, 'sig_corr_unconnected':sig_corr_unconnected, 'normal_receptor_list': normal_receptor_list, 'cancer_ligand_list': cancer_ligand_list,'mean_normal_receptor':mean_normal_receptor, 'mean_cancer_ligand':mean_cancer_ligand,'normal_receptor_list_uc': normal_receptor_list_uc, 'cancer_ligand_list_uc': cancer_ligand_list_uc, 'n_con_cancer':n_con_cancer, 'n_con_normal': n_con_normal}

def lig_rec_plot_normal_cancer_NC(adata,cell_subtype,impute=True,lig_rec_pair=('PDGFB', 'LRP1'),all_genes=False,corr_coef=spearmanr):
    if (impute):
        #print('Using imputed values')
        adata.X=adata.obsm['X_normalized_scVI']
        assert(np.sum(np.asarray(adata.X))==np.sum(np.asarray(adata.obsm['X_normalized_scVI'])))
    else:
        print('Not using imputed values')
    print('Considering '+Normal_cells[cell_subtype])
    print('Cancer_cell')
    print(Cancer_cell)
    adata_cancer_subtype=adata[(adata.obs[ct_level]==(Cancer_cell)) | (adata.obs[ct_level]==(Normal_cells[cell_subtype]))]
    if 'connectivities' in adata_cancer_subtype.obsp:
        connectivities_empty = adata_cancer_subtype.obsp['connectivities'].nnz == 0  # Check if the number of non-zero elements is 0
    else:
        connectivities_empty = True
    # If 'connectivities' is empty, set it to 'spatial_connectivities'
    if connectivities_empty:
        print('Note that connectivities was empty, setting it to spatial_connectivities')
        adata_cancer_subtype.obsp['connectivities'] = adata_cancer_subtype.obsp['spatial_connectivities']
    assert(len(adata_cancer_subtype)==adata_cancer_subtype.obsp['connectivities'].shape[0])
    assert(len(adata_cancer_subtype)==adata_cancer_subtype.obsp['connectivities'].shape[1])
    adjacency_matrix=adata_cancer_subtype.obsp['connectivities'].toarray()
    adjacency_matrix_df=pd.DataFrame(adjacency_matrix,index=adata_cancer_subtype.obs.index, columns=adata_cancer_subtype.obs.index)
    index_list1 = adata_cancer_subtype.obs[(adata_cancer_subtype.obs[ct_level] == Cancer_cell)].index
    index_list2 = adata_cancer_subtype.obs[(adata_cancer_subtype.obs[ct_level]==Normal_cells[cell_subtype])].index
    index_list1 = index_list1.astype(str)
    index_list2 = index_list2.astype(str)
    adjacency_matrix_rows = adjacency_matrix_df.loc[index_list2, index_list1]
    row_labels = adjacency_matrix_rows.index[np.where(adjacency_matrix_rows == 1.0)[0]]
    col_labels = adjacency_matrix_rows.columns[np.where(adjacency_matrix_rows == 1.0)[1]]
    row_labels_uc = adjacency_matrix_rows.index[np.where((adjacency_matrix_rows != 1.0).all(axis=1) & (adjacency_matrix_rows.sum(axis=1) == 0))[0]]
    col_labels_uc = adjacency_matrix_rows.columns[np.where((adjacency_matrix_rows != 1.0).all(axis=0) & (adjacency_matrix_rows.sum(axis=0) == 0))[0]]
    ### For NC direction
    normal_ligand=lig_rec_pair[0]
    cancer_receptor=lig_rec_pair[1]
    idx_ligand=adata.var_names.str.upper().get_loc(normal_ligand)
    idx_receptor=adata.var_names.str.upper().get_loc(cancer_receptor)
    normal_ligand_list = adata_cancer_subtype[row_labels, :].X[:, idx_ligand]
    cancer_receptor_list = adata_cancer_subtype[col_labels, :].X[:, idx_receptor]
    mean_cancer_list=np.mean(cancer_receptor_list)
    mean_normal_list=np.mean(normal_ligand_list)
    corr_connected=round(corr_coef(normal_ligand_list, cancer_receptor_list)[0], 4)
    print('corr_connected NC:',corr_connected)
    sig_corr_connected=round(corr_coef(normal_ligand_list, cancer_receptor_list)[1], 4)
    connected_df_NC=pd.DataFrame()
    connected_df_NC['normal_ligand_list']=normal_ligand_list
    connected_df_NC['cancer_receptor_list']=cancer_receptor_list
    n_con_cancer=len(cancer_receptor_list)
    n_con_normal=len(normal_ligand_list)
    mean_cancer_receptor=np.mean(cancer_receptor_list)
    mean_normal_ligand=np.mean(normal_ligand_list)
    max_cancer_receptor=np.max(cancer_receptor_list)
    max_normal_ligand=np.max(normal_ligand_list)
    min_cancer_receptor=np.min(cancer_receptor_list)
    min_normal_ligand=np.min(normal_ligand_list)
    normal_ligand_list_uc = adata_cancer_subtype[row_labels_uc,:].X[:,idx_ligand]
    cancer_receptor_list_uc = adata_cancer_subtype[col_labels_uc,:].X[:,idx_receptor]
    normal_ligand_list_uc_repeats=np.repeat(normal_ligand_list_uc,len(cancer_receptor_list_uc))
    cancer_receptor_list_uc_repeats=np.repeat(cancer_receptor_list_uc,len(normal_ligand_list_uc))
    cancer_receptor_list_uc = np.array(cancer_receptor_list_uc)
    normal_ligand_list_uc = np.array(normal_ligand_list_uc)
    N = 0.1*len(normal_ligand_list_uc)*len(cancer_receptor_list_uc)
    sampled_indices1 = np.random.choice(len(normal_ligand_list_uc), int(N))
    sampled_indices2 = np.random.choice(len(cancer_receptor_list_uc), int(N))
    sampled_normal_ligand_list_uc = normal_ligand_list_uc[sampled_indices1]
    sampled_cancer_receptor_list_uc = cancer_receptor_list_uc[sampled_indices2]
    unconnected_df_NC=pd.DataFrame()
    unconnected_df_NC['cancer_receptor_list_uc']=sampled_cancer_receptor_list_uc
    unconnected_df_NC['normal_ligand_list_uc']=sampled_normal_ligand_list_uc
    n_uncon_cancer=len(sampled_cancer_receptor_list_uc)
    n_uncon_normal=len(sampled_normal_ligand_list_uc)
    corr_unconnected=round(corr_coef(sampled_cancer_receptor_list_uc, sampled_normal_ligand_list_uc)[0], 4)
    sig_corr_unconnected=round(corr_coef(sampled_cancer_receptor_list_uc, sampled_normal_ligand_list_uc)[1], 4)
    if ((corr_connected > 0.05) or (corr_connected < -0.05)):
        print('Plotting scatterplot for potential ligand-receptor pair')
        plt.close()
        plt.figure(figsize=(10, 5))
        plt.subplot(1, 2, 1)
        plt.violinplot(normal_ligand_list,showmedians=True)
            #plt.text(0.5, -0.4, lig_rec_pair[0], horizontalalignment='center', verticalalignment='center', fontsize=12)
        plt.text(0.5, min(normal_ligand_list) - 0.1, lig_rec_pair[0], horizontalalignment='center', verticalalignment='center', fontsize=12)
        plt.subplot(1, 2, 2)
        plt.violinplot(cancer_receptor_list,showmedians=True)
        plt.text(0.5, min(cancer_receptor_list) - 0.1, lig_rec_pair[1], horizontalalignment='center', verticalalignment='center', fontsize=12)
        plt.savefig(path+Normal_cells[cell_subtype]+'_imgs/'+'Violin_plot_'+Normal_cells[cell_subtype]+'_'+ct_level+ '_'+ version+'_'+lig_rec_pair[0]+'-'+lig_rec_pair[1]+'.png')
        plt.close()
        plt.hist2d(normal_ligand_list, cancer_receptor_list, bins=30, cmap='viridis')
        plt.colorbar(label='Frequency')
        plt.xlabel(lig_rec_pair[0])
        plt.ylabel(lig_rec_pair[1])
        plt.title('2D Histogram of ligand-receptor expression')
        plt.savefig(path+Normal_cells[cell_subtype]+'_imgs/'+'2D_hist_'+Normal_cells[cell_subtype]+'_'+ct_level+ '_'+ version+'_'+lig_rec_pair[0]+'-'+lig_rec_pair[1]+'.png')
        plt.close()
        plt.figure(figsize=(10, 5))
        try:
            sns.kdeplot(data=connected_df_NC, color="blue")
        except ValueError as e:
            print("Error occurred while plotting KDE:", e)
            print("Gene pair responsible is:\n")
            print(lig_rec_pair[0]+' - '+lig_rec_pair[1])
        #plt.close()
        plt.title(Normal_cells[cell_subtype]+'->C'+' '+lig_rec_pair[0]+' - '+lig_rec_pair[1])
        plt.scatter(normal_ligand_list,cancer_receptor_list,color='blue',label='Connected pairs, correlation= '+str(corr_connected), alpha=0.5)
        plt.xlabel(Normal_cells[cell_subtype]+ ' mean = '+str(np.mean(normal_ligand_list)))
        plt.ylabel('Epithelial '+' mean = '+str(np.mean(cancer_receptor_list)))
        plt.legend()
        if (all_genes):
            print('Doing all genes!')
            if ((corr_connected>0.8) or (corr_connected < -0.8)):
                plt.savefig(path+Normal_cells[cell_subtype]+'_imgs/'+'NC_Connected_vs_unconnected_pair_'+Normal_cells[cell_subtype]+'_'+ct_level+ '_'+ version+'_'+lig_rec_pair[0]+'-'+lig_rec_pair[1]+'_normal_ligand.png')
        else:
            plt.savefig(path+Normal_cells[cell_subtype]+'_imgs/'+'NC_Connected_vs_unconnected_pair_'+Normal_cells[cell_subtype]+'_'+ct_level+ '_'+ version+'_'+lig_rec_pair[0]+'-'+lig_rec_pair[1]+'_normal_ligand.png')
    return {'corr_connected':corr_connected, 'sig_corr_connected':sig_corr_connected, 'corr_unconnected':corr_unconnected, 'sig_corr_unconnected':sig_corr_unconnected, 'normal_ligand_list': normal_ligand_list, 'cancer_receptor_list': cancer_receptor_list,'mean_cancer_receptor':mean_cancer_receptor, 'mean_normal_ligand':mean_normal_ligand,'normal_ligand_list_uc': normal_ligand_list_uc, 'cancer_receptor_list_uc': cancer_receptor_list_uc, 'n_con_cancer':n_con_cancer, 'n_con_normal': n_con_normal}


h5ad_file='anndata_GR_with_delaunay_40_microns_'+ct_level+'_'+version+'.h5ad'

adata=ad.read_h5ad(path+h5ad_file)

print('Read in the following anndata object '+ str(adata))

if (args.ct_level=='ct_level1') or (args.ct_level=='ct_level0'): #or ct_level=='annot_level0'):
    Cancer_cell='Epithelial'
    print('Doing Epithelial cell')
if args.ct_level=='annot_level0':
    Cancer_cell='Epithelial cells'
    print('Doing Epithelial cells')
if args.ct_level=='ct_clean':
    Cancer_cell=='Cancer_cell'
    print('Doing Cancer_cell, not Epithelial cell')

print('Cancer cell type is '+Cancer_cell)

Normal_cells = np.unique(adata.obs[ct_level])[[x!=Cancer_cell for x in np.unique(adata.obs[ct_level])]]

print('Normal cells:\n')
print(Normal_cells)

cell_idx=np.where(Normal_cells==cell_subtype)[0][0]

print('cell_idx')
print(cell_idx)

directory_name=path+Normal_cells[cell_idx]+'_imgs/'

print('directory_name')
print(directory_name)

# Check if the directory exists
if not os.path.exists(directory_name):
    # Create the directory if it doesn't exist
    os.makedirs(directory_name)

ligand_receptor_df=pd.read_csv("/path/to/file/Human-2020-Jin-LR-pairs.csv")

ligand_receptor_df=ligand_receptor_df[(ligand_receptor_df['annotation']=='Cell-Cell Contact') | (ligand_receptor_df['annotation']=='ECM-Receptor')]

# Split ligands and receptors if they contain entries separated by '_'
ligand_receptor_df[['ligand1', 'ligand2']] = ligand_receptor_df['ligand'].str.split('_', n=1, expand=True)
ligand_receptor_df[['receptor1', 'receptor2']] = ligand_receptor_df['receptor'].str.split('_', n=1, expand=True)

ligand_receptor_df_subset=ligand_receptor_df[['ligand1','ligand2','receptor1','receptor2','pathway_name']]

ligand_receptor_df_subset = ligand_receptor_df_subset.applymap(lambda x: x.upper() if isinstance(x, str) else x)

panel_genes = pd.Series(adata.var_names)

panel_genes = panel_genes.str.upper()

if (all_genes):
    print('Generating pairwise combinations for all genes')
    # Filter out genes with the name "Blank" or "BLANK"
    filtered_genes = panel_genes[~panel_genes.isin(['Blank', 'BLANK'])]

    # Generate all pairwise combinations
    gene_pairs = list(itertools.product(filtered_genes, repeat=2))

    # Define the path for the output file
    output_file_path = path + "gene_pairs_all.txt"

    # Write pairwise combinations to a file
    with open(output_file_path, "w") as file:
        for pair in gene_pairs:
            file.write(" ".join(pair) + "\n")

    print(f'Pairwise combinations of genes written to {output_file_path}')

# Initialize an empty list to store valid pairs
valid_pairs = []
pathway_list = []

# Iterate through rows
for index, row in ligand_receptor_df_subset.iterrows():
    # Iterate through combinations
    for ligand, receptor, pathway in [(row['ligand1'], row['receptor1'], row['pathway_name']),
                             (row['ligand1'], row['receptor2'], row['pathway_name']),
                             (row['ligand2'], row['receptor1'], row['pathway_name']),
                             (row['ligand2'], row['receptor2'], row['pathway_name'])]:
        # Check if both ligand and receptor are in panel_genes
        if ligand in list(panel_genes) and receptor in list(panel_genes):
            valid_pairs.append((ligand, receptor, pathway))
            #pathway_list.append(pathway)

print('valid_pairs:\n')
print(valid_pairs)

with open(path+"gene_pairs_known_ligrec.txt", "w") as file:
    for pair in valid_pairs:
        file.write(" ".join(pair) + "\n")
        
print("Wrote down known ligrec pairs with pathway")


if (all_genes):
    print('Note: reading in all gene pairs')
    #file_path='/path/to/file/all_gene_pairs_TNBC.txt'
    file_path=path+'gene_pairs_all.txt'
if (ligrec_known):
    print('Note:reading in known ligrec pairs')
    file_path=path+'gene_pairs_known_ligrec.txt'
    #pathways_path=path+'pathways_for_gene_pairs.txt'
if (DE_known):
    print('Note: reading in DE gene pairs')
    file_path=path+cell_subtype+'_gene_pairs_'+'_'+ct_level+'_'+version+'_'+sample+'_cat.txt'

print('The gene combinations file is:\n')
 
print(file_path)
 
if (all_genes or DE_known):
    print('Doing all_genes or DE_known')
    file_contents=[]
    with open(file_path, 'r') as file:
    # Read each line in the file
        lines = file.readlines()
        for line in lines:
            genes = line.strip().split()
            print('Reading in all gene-gene pairs or DE genes')
            gene1, gene2 = genes[0],genes[1]
            file_contents.append((gene1, gene2))
        print('file_contents:\n')
        print(file_contents)

if (all_genes):
   print('Looking at all pairwise combinations')
   valid_pairs=file_contents
   #print('Length:\n')
   #print(len(valid_pairs))
if (ligrec_known):
   print('Looking at known ligand-receptor pairs')
   valid_pairs=valid_pairs
if (DE_known):
   print('Looking at DE pairs')

print('Length of gene-gene pairs:\n')
print(len(valid_pairs))
 
combination=valid_pairs

print('Length of split combination:\n')
print(len(combination))

print('start_index')
print(start_index)
print('end_index')
print(end_index)

def process_split_known_ligrec(combination,start,end):
    data_CN=[]
    data_NC=[]
    pathway_list=[]
    #if all_genes:
    #    print('Looking at all ligand-receptor pairs')
    #if ligrec_known:
    print('Looking at known ligand-receptor pairs')
    pathway_list=[]
    #if DE_known:
    #print('Looking at ligand-receptor pairs from DE analysis')
    for j in range(start,end):
        print(f"gene pair {combination[j][0]} {combination[j][1]}")
        print(f"pathway {combination[j][2]}")
        pathway_list.append(combination[j][2])
        #corr_connected_CN, sig_corr_connected_CN, n_con_cancer, n_con_normal, mean_cancer_ligand, mean_normal_receptor, max_cancer_ligand, max_normal_receptor, min_cancer_ligand, min_normal_receptor=lig_rec_plot_normal_cancer_CN(adata,cell_subtype=cell_idx,impute=True,lig_rec_pair=combination[j],all_genes=all_genes,corr_coef=correlation_function)
        #corr_connected_NC, sig_corr_connected_NC, n_con_cancer, n_con_normal, mean_normal_ligand, mean_cancer_receptor, max_normal_ligand, max_cancer_receptor, min_normal_ligand, min_cancer_receptor=lig_rec_plot_normal_cancer_NC(adata,cell_subtype=cell_idx,impute=True,lig_rec_pair=combination[j],all_genes=all_genes,corr_coef=correlation_function)
        results_CN=lig_rec_plot_normal_cancer_CN(adata,cell_subtype=cell_idx,impute=True,lig_rec_pair=combination[j],all_genes=all_genes,corr_coef=correlation_function)
        results_NC=lig_rec_plot_normal_cancer_NC(adata,cell_subtype=cell_idx,impute=True,lig_rec_pair=combination[j],all_genes=all_genes,corr_coef=correlation_function)
        data_CN.append([combination[j][0], combination[j][1],combination[j][2],results_CN['corr_connected'],results_CN['sig_corr_connected'],results_CN['corr_unconnected'],results_CN['sig_corr_unconnected'],
        results_CN['normal_receptor_list'].tolist(),results_CN['cancer_ligand_list'].tolist(),results_CN['mean_normal_receptor'],results_CN['mean_cancer_ligand'],results_CN['normal_receptor_list_uc'],results_CN['cancer_ligand_list_uc'],results_CN['n_con_cancer'],results_CN['n_con_normal']])
        data_NC.append([combination[j][0], combination[j][1],combination[j][2],results_NC['corr_connected'],results_NC['sig_corr_connected'],results_NC['corr_unconnected'],results_NC['sig_corr_unconnected'],
        results_NC['normal_ligand_list'].tolist(),results_NC['cancer_receptor_list'].tolist(),results_NC['mean_cancer_receptor'],results_NC['mean_normal_ligand'],results_NC['normal_ligand_list_uc'],results_NC['cancer_receptor_list_uc'],results_NC['n_con_cancer'],results_NC['n_con_normal']])
    return data_CN, data_NC

#corr_connected, sig_corr_connected,corr_unconnected, sig_corr_unconnected, 
#normal_receptor_list, cancer_ligand_list,mean_normal_receptor, mean_cancer_ligand, 
#normal_receptor_list_uc, cancer_ligand_list_uc,n_con_cancer, n_con_normal

#normal_ligand_list_uc, cancer_receptor_list_uc

def process_split(combination,start,end):
    data_CN=[]
    data_NC=[]
    print('Looking at all genes')
    for j in range(start,end):
        print(f"gene pair {combination[j][0]} {combination[j][1]}")
        #corr_connected_CN, sig_corr_connected_CN, n_con_cancer, n_con_normal, mean_cancer_ligand, mean_normal_receptor, max_cancer_ligand, max_normal_receptor, min_cancer_ligand, min_normal_receptor=lig_rec_plot_normal_cancer_CN(adata,cell_subtype=cell_idx,impute=True,lig_rec_pair=combination[j],all_genes=all_genes,corr_coef=correlation_function)
        #corr_connected_NC, sig_corr_connected_NC, n_con_cancer, n_con_normal, mean_normal_ligand, mean_cancer_receptor, max_normal_ligand, max_cancer_receptor, min_normal_ligand, min_cancer_receptor=lig_rec_plot_normal_cancer_NC(adata,cell_subtype=cell_idx,impute=True,lig_rec_pair=combination[j],all_genes=all_genes,corr_coef=correlation_function)
        results_CN=lig_rec_plot_normal_cancer_CN(adata,cell_subtype=cell_idx,impute=True,lig_rec_pair=combination[j],all_genes=all_genes,corr_coef=correlation_function)
        results_NC=lig_rec_plot_normal_cancer_NC(adata,cell_subtype=cell_idx,impute=True,lig_rec_pair=combination[j],all_genes=all_genes,corr_coef=correlation_function)
        data_CN.append([combination[j][0], combination[j][1],results_CN['corr_connected'],results_CN['sig_corr_connected'],results_CN['corr_unconnected'],results_CN['sig_corr_unconnected'],
        results_CN['normal_receptor_list'].tolist(),results_CN['cancer_ligand_list'].tolist(),results_CN['mean_normal_receptor'],results_CN['mean_cancer_ligand'],results_CN['normal_receptor_list_uc'],results_CN['cancer_ligand_list_uc'],
        results_CN['normal_receptor_list_uc'],results_CN['cancer_ligand_list_uc'],results_CN['n_con_cancer'],results_CN['n_con_normal']])
        data_NC.append([combination[j][0], combination[j][1],results_NC['corr_connected'],results_NC['sig_corr_connected'],results_NC['corr_unconnected'],results_NC['sig_corr_unconnected'],
        results_NC['normal_ligand_list'].tolist(),results_NC['cancer_receptor_list'].tolist(),results_NC['mean_cancer_receptor'],results_NC['mean_normal_ligand'],
        results_NC['normal_ligand_list_uc'],results_NC['cancer_receptor_list_uc'],results_NC['n_con_cancer'],results_NC['n_con_normal']])
    return data_CN, data_NC


# Process the split
# Process the split
if (ligrec_known):
     data_CN1, data_NC1 =process_split_known_ligrec(combination, start_index, end_index)
     df_CN = pd.DataFrame(data_CN1, columns=["Ligand", "Receptor", "Pathway", "corr_connected", "sig_corr_connected", "corr_unconnected","sig_corr_unconnected","normal_receptor_list","cancer_ligand_list","mean_normal_receptor","mean_cancer_ligand", "normal_receptor_list_uc","cancer_ligand_list_uc","n_con_cancer","n_con_normal"])
     df_NC = pd.DataFrame(data_NC1, columns=["Ligand", "Receptor", "Pathway", "corr_connected", "sig_corr_connected", "corr_unconnected","sig_corr_unconnected","normal_ligand_list","cancer_receptor_list","mean_cancer_receptor","mean_normal_ligand", "normal_ligand_list_uc","cancer_receptor_list_uc","n_con_cancer","n_con_normal"])
if (all_genes):
     data_CN1, data_NC1 =process_split(combination, start_index, end_index)
     df_CN = pd.DataFrame(data_CN1, columns=["Ligand", "Receptor", "corr_connected", "sig_corr_connected", "corr_unconnected","sig_corr_unconnected","normal_receptor_list","cancer_ligand_list","mean_normal_receptor","mean_cancer_ligand", "normal_receptor_list_uc","cancer_ligand_list_uc","n_con_cancer","n_con_normal"])
     df_NC = pd.DataFrame(data_NC1, columns=["Ligand", "Receptor", "corr_connected", "sig_corr_connected", "corr_unconnected","sig_corr_unconnected","normal_ligand_list","cancer_receptor_list","mean_cancer_receptor","mean_normal_ligand", "normal_ligand_list_uc","cancer_receptor_list_uc","n_con_cancer","n_con_normal"])

print('all_genes')
print(all_genes)

#df_CN = pd.DataFrame(data_CN1, columns=["Ligand", "Receptor","Pathway", "corr_connected", "sig_corr_connected", "corr_unconnected","sig_corr_unconnected","normal_receptor_list","cancer_ligand_list","mean_normal_receptor","mean_cancer_ligand", "n_con_cancer","n_con_normal"])

#df_NC = pd.DataFrame(data_NC1, columns=["Ligand", "Receptor","Pathway", "corr_connected", "sig_corr_connected", "corr_unconnected","sig_corr_unconnected","normal_ligand_list","cancer_receptor_list","mean_cancer_receptor","mean_normal_ligand", "n_con_cancer","n_con_normal"])

print('df_CN.shape')
print(df_CN.shape)

print('df_NC.shape')
print(df_NC.shape)

#print('pathway_list')
#print(pathway_list)

#print('len(pathway_list)')
#print(len(pathway_list))

#df_CN['pathway']=pathway_list
#df_NC['pathway']=pathway_list

print(df_CN)

print(df_NC)
        
df_CN.to_csv(path+Normal_cells[cell_idx]+'_imgs/'+Normal_cells[cell_idx]+'_'+ct_level+ '_'+ version+'_CN_'+'_main'+'_'+all_flag+'_'+'part_'+str(split_part)+'_complete_connected.csv', index=False)
df_NC.to_csv(path+Normal_cells[cell_idx]+'_imgs/'+Normal_cells[cell_idx]+'_'+ct_level+ '_'+ version+'_NC_'+'_main'+'_'+all_flag+'_'+'part_'+str(split_part)+'_complete_connected.csv', index=False)

print('Finished run successfully')

# Record end time
end_time = time.time()

# Calculate elapsed time
elapsed_time = end_time - start_time
print(f"Elapsed time: {elapsed_time} seconds")




