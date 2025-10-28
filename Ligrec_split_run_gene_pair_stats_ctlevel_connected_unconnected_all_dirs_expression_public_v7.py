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
import scipy.sparse as sp

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
parser.add_argument('--hfad_file', type=str, required=True, help='h5ad file')
parser.add_argument('--subset', action='store_true', help='Consider a subset of the original data')
parser.add_argument('--subset_type', type=str, choices=['60p', '70p', '80p'], help='Which subset to use: 60p, 70p, or 80p')

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
if args.subset:
    print("Considering data subset")
    subset=True

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
hfad_file=args.hfad_file
subset=args.subset
subset_type = args.subset_type

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
print('hfad file:\n')
print(hfad_file)
print('subset:\n')
print(subset)
print('subset_type:\n')
print(subset_type)

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

# Determine which AnnData file to load based on the subset argument
if subset:
    print("Considering data subset")
    
    # Ensure that the subset_type is provided
    if not subset_type:
        raise ValueError("Please provide the subset type (60p, 70p, or 80p) using the --subset_type argument.")
    
    h5ad_file = f'anndata_subset_{ct_level}_{version}_{subset_type}.h5ad'
else:
    h5ad_file = hfad_file

# Load the appropriate AnnData object
adata = ad.read_h5ad(os.path.join(path, h5ad_file))
print(f'Read in the following anndata object: {str(adata)}')

if isinstance(adata.obs.index, pd.CategoricalIndex):
    adata.obs.index = adata.obs.index.astype(str)

for col in adata.obs.select_dtypes(['category']).columns:
    adata.obs[col] = adata.obs[col].astype(str)


def lig_rec_plot_normal_cancer_CN(adata,cell_subtype,impute=True,lig_rec_pair=('PDGFB', 'LRP1'),all_genes=False):
    if (impute):
        #print('Using imputed values')
        adata.X=adata.obsm['X_normalized_scVI']
        #adata.X=adata.obsm['X_scVI']
        #assert(np.sum(np.asarray(adata.X))==np.sum(np.asarray(adata.obsm['X_normalized_scVI'])))
    else:
        print('Not using imputed values')
    print('Considering '+Normal_cells[cell_subtype])
    print('Cancer_cell')
    print(Cancer_cell)
    adata_cancer_subtype=adata[(adata.obs[ct_level]==(Cancer_cell)) | (adata.obs[ct_level]==(Normal_cells[cell_subtype]))]
    #print('adata_cancer_subtype.obsp[spatial_connectivities]:\n')
    #print(adata_cancer_subtype.obsp['spatial_connectivities'])
    # Check if 'connectivities' exists in obsp and if it's empty
    if 'connectivities' in adata_cancer_subtype.obsp:
        connectivities_empty = adata_cancer_subtype.obsp['connectivities'].nnz == 0  # Check if the number of non-zero elements is 0
    else:
        connectivities_empty = True
    # If 'connectivities' is empty, set it to 'spatial_connectivities'
    if connectivities_empty:
        #print('Note that connectivities was empty, setting it to spatial_connectivities')
        adata_cancer_subtype.obsp['connectivities'] = adata_cancer_subtype.obsp['spatial_connectivities']
    assert(len(adata_cancer_subtype)==adata_cancer_subtype.obsp['connectivities'].shape[0])
    assert(len(adata_cancer_subtype)==adata_cancer_subtype.obsp['connectivities'].shape[1])
    adjacency_matrix=adata_cancer_subtype.obsp['connectivities'].toarray()
    #adjacency_matrix=adata_cancer_subtype.obsp['connectivities']
    #adjacency_matrix_df = pd.DataFrame.sparse.from_spmatrix(adjacency_matrix, index=adata_cancer_subtype.obs.index, columns=adata_cancer_subtype.obs.index)
    adjacency_matrix_df=pd.DataFrame(adjacency_matrix,index=adata_cancer_subtype.obs.index, columns=adata_cancer_subtype.obs.index)
    index_list1 = adata_cancer_subtype.obs[(adata_cancer_subtype.obs[ct_level] == Cancer_cell)].index
    index_list2 = adata_cancer_subtype.obs[(adata_cancer_subtype.obs[ct_level]==Normal_cells[cell_subtype])].index
    index_list1 = index_list1.astype(str)
    index_list2 = index_list2.astype(str)
    adjacency_matrix_rows = adjacency_matrix_df.loc[index_list2, index_list1]
    #print('adjacency_matrix_rows:\n')
    #print(adjacency_matrix_rows)
    #print('adjacency_matrix_rows.shape')
    #print(adjacency_matrix_rows.shape)
    row_labels = adjacency_matrix_rows.index[np.where(adjacency_matrix_rows == 1.0)[0]]
    col_labels = adjacency_matrix_rows.columns[np.where(adjacency_matrix_rows == 1.0)[1]]
    row_labels_uc = adjacency_matrix_rows.index[np.where((adjacency_matrix_rows != 1.0).all(axis=1) & (adjacency_matrix_rows.sum(axis=1) == 0))[0]]
    col_labels_uc = adjacency_matrix_rows.columns[np.where((adjacency_matrix_rows != 1.0).all(axis=0) & (adjacency_matrix_rows.sum(axis=0) == 0))[0]]
    #row_labels = adjacency_matrix_rows.index[adjacency_matrix_rows.sparse.to_dense() == 1.0]
    #col_labels = adjacency_matrix_rows.columns[adjacency_matrix_rows.sparse.to_dense() == 1.0]
    #row_labels_uc = adjacency_matrix_rows.index[(adjacency_matrix_rows.sparse.to_dense() != 1.0).all(axis=1) & (adjacency_matrix_rows.sparse.to_dense().sum(axis=1) == 0)]
    #col_labels_uc = adjacency_matrix_rows.columns[(adjacency_matrix_rows.sparse.to_dense() != 1.0).all(axis=0) & (adjacency_matrix_rows.sparse.to_dense().sum(axis=0) == 0)]
    ### For CN direction
    cancer_ligand=lig_rec_pair[0]
    normal_receptor=lig_rec_pair[1]
    try:
        idx_ligand=adata.var_names.str.upper().get_loc(cancer_ligand)
        idx_receptor=adata.var_names.str.upper().get_loc(normal_receptor)
        cancer_ligand_list = adata_cancer_subtype[row_labels, :].X[:, idx_ligand]
        normal_receptor_list = adata_cancer_subtype[col_labels, :].X[:, idx_receptor]
        n_con_cancer=len(cancer_ligand_list)
        n_con_normal=len(normal_receptor_list)
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
    except KeyError:
        print(cancer_ligand + ' and ' + normal_receptor + ' are not found in the data, trying the next one')
        return {
            'normal_receptor_list': np.array([np.nan]),
            'cancer_ligand_list': np.array([np.nan]),
            'normal_receptor_list_uc': np.array([np.nan]),
            'cancer_ligand_list_uc': np.array([np.nan]),
            'n_con_cancer': np.nan,
            'n_con_normal': np.nan
        }
    #connected_df_CN=pd.DataFrame()
    #connected_df_CN['normal_receptor_list']=normal_receptor_list
    #connected_df_CN['cancer_ligand_list']=cancer_ligand_list
    return {'normal_receptor_list': normal_receptor_list, 'cancer_ligand_list': cancer_ligand_list,'normal_receptor_list_uc': normal_receptor_list_uc, 'cancer_ligand_list_uc': cancer_ligand_list_uc, 'n_con_cancer':n_con_cancer, 'n_con_normal': n_con_normal}

#, 'mean_cancer_ligand':mean_cancer_ligand, 'mean_normal_receptor':mean_normal_receptor

def lig_rec_plot_normal_cancer_NC(adata,cell_subtype,impute=True,lig_rec_pair=('PDGFB', 'LRP1'),all_genes=False):
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
        #print('Note that connectivities was empty, setting it to spatial_connectivities')
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
    try:
        idx_ligand=adata.var_names.str.upper().get_loc(normal_ligand)
        idx_receptor=adata.var_names.str.upper().get_loc(cancer_receptor)
        normal_ligand_list = adata_cancer_subtype[row_labels, :].X[:, idx_ligand]
        cancer_receptor_list = adata_cancer_subtype[col_labels, :].X[:, idx_receptor]
        n_con_cancer=len(cancer_receptor_list)
        n_con_normal=len(normal_ligand_list)
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
    except KeyError:
        print(normal_ligand + ' and ' + cancer_receptor + ' are not found in the data, trying the next one')
        #continue
        return {
            'normal_ligand_list': np.array([np.nan]),
            'cancer_receptor_list': np.array([np.nan]),
            'normal_ligand_list_uc': np.array([np.nan]),
            'cancer_receptor_list_uc': np.array([np.nan]),
            'n_con_cancer': np.nan,
            'n_con_normal': np.nan
        }
    #corr_unconnected=round(corr_coef(sampled_cancer_receptor_list_uc, sampled_normal_ligand_list_uc)[0], 4)
    #sig_corr_unconnected=round(corr_coef(sampled_cancer_receptor_list_uc, sampled_normal_ligand_list_uc)[1], 4)
    return {'normal_ligand_list': normal_ligand_list, 'cancer_receptor_list': cancer_receptor_list,'normal_ligand_list_uc': normal_ligand_list_uc, 'cancer_receptor_list_uc': cancer_receptor_list_uc, 'n_con_cancer':n_con_cancer, 'n_con_normal': n_con_normal}

def split_anndata_by_coordinates_into_9(adata, x_key='x', y_key='y'):
    # Get the x and y coordinates
    x_coords = adata.obs[x_key].values
    y_coords = adata.obs[y_key].values

    # Calculate the 1/3 and 2/3 quantiles to split the data
    x_split1 = np.quantile(x_coords, 1/3)
    x_split2 = np.quantile(x_coords, 2/3)
    y_split1 = np.quantile(y_coords, 1/3)
    y_split2 = np.quantile(y_coords, 2/3)

    # Create masks for each region
    masks = {
        'region_1': (x_coords <= x_split1) & (y_coords > y_split2),
        'region_2': (x_coords > x_split1) & (x_coords <= x_split2) & (y_coords > y_split2),
        'region_3': (x_coords > x_split2) & (y_coords > y_split2),
        'region_4': (x_coords <= x_split1) & (y_coords > y_split1) & (y_coords <= y_split2),
        'region_5': (x_coords > x_split1) & (x_coords <= x_split2) & (y_coords > y_split1) & (y_coords <= y_split2),
        'region_6': (x_coords > x_split2) & (y_coords > y_split1) & (y_coords <= y_split2),
        'region_7': (x_coords <= x_split1) & (y_coords <= y_split1),
        'region_8': (x_coords > x_split1) & (x_coords <= x_split2) & (y_coords <= y_split1),
        'region_9': (x_coords > x_split2) & (y_coords <= y_split1)
    }

    # Slice the AnnData object into 8 parts
    adata_parts = {key: adata[mask].copy() for key, mask in masks.items()}

    return adata_parts
    
# Example usage
# Assume `adata` is your AnnData object and `obs` contains 'x' and 'y' coordinates
adata_parts = split_anndata_by_coordinates_into_9(adata,x_key='center_x',y_key='center_y')

# Now you have adata_q1, adata_q2, adata_q3, and adata_q4 as separate AnnData objects

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

if (all_genes):
    print('Note: reading in all gene pairs')
    file_path='/path/to/file/all_gene_pairs_TNBC.txt'
    #file_path=path+'all_gene_pairs_'+ct_level+'_'+version+'_'+sample+'.txt'
if (ligrec_known):
    print('Note: reading in known ligrec pairs')
    file_path=path+'gene_pairs_known_ligrec_Ovarian.txt'
    #pathways_path=path+'pathways_for_gene_pairs.txt'
if (DE_known):
    print('Note: reading in DE gene pairs')
    file_path=path+cell_subtype+'_gene_pairs_'+'_'+ct_level+'_'+version+'_'+sample+'_cat.txt'

print('file path is:\n')
print(file_path)
#print('The gene combinations file is:\n')
 
#print(file_path)
 
file_contents=[]
if (ligrec_known):
    with open(file_path, 'r') as file:
        print('Reading in known lig-rec pairs with pathway')
    # Read each line in the file
        lines = file.readlines()
        for line in lines:
            genes = line.strip().split()
            gene1, gene2, pathway = genes[0],genes[1], genes[2]
            file_contents.append((gene1, gene2, pathway))
        print('file_contents:\n')
        print(file_contents)

if (all_genes):
   print('Looking at all pairwise combinations')
   valid_pairs=file_contents
   #print('Length:\n')
   #print(len(valid_pairs))
if (ligrec_known):
   print('Looking at known ligand-receptor pairs')
   valid_pairs=file_contents
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

def process_split_known_ligrec(combination, start, end, adata_parts):
    data_CN = []
    data_NC = []
    pathway_list = []

    print('Looking at known ligand-receptor pairs')
    for j in range(start, end):
        print(f"gene pair {combination[j][0]} {combination[j][1]}")
        print(f"pathway {combination[j][2]}")
        pathway_list.append(combination[j][2])

        # Iterate over all regions
        for region_key, adata_region in adata_parts.items():
            print(f"Processing {region_key}")

            results_CN = lig_rec_plot_normal_cancer_CN(
                adata_region,
                cell_subtype=cell_idx,
                impute=True,
                lig_rec_pair=combination[j],
                all_genes=all_genes,
            )
            results_NC = lig_rec_plot_normal_cancer_NC(
                adata_region,
                cell_subtype=cell_idx,
                impute=True,
                lig_rec_pair=combination[j],
                all_genes=all_genes
            )

            data_CN.append([
                combination[j][0], combination[j][1], combination[j][2],
                results_CN['normal_receptor_list'].tolist(), results_CN['cancer_ligand_list'].tolist(),
                results_CN['normal_receptor_list_uc'], results_CN['cancer_ligand_list_uc'],
                #results_CN['mean_cancer_ligand'],results_CN['mean_normal_receptor'],
                results_CN['n_con_cancer'], results_CN['n_con_normal']
            ])
            data_NC.append([
                combination[j][0], combination[j][1], combination[j][2],
                results_NC['normal_ligand_list'].tolist(), results_NC['cancer_receptor_list'].tolist(),
                results_NC['normal_ligand_list_uc'], results_NC['cancer_receptor_list_uc'],
                #results_NC['mean_normal_ligand'],results_NC['mean_cancer_receptor'],
                results_NC['n_con_cancer'], results_NC['n_con_normal']
            ])
    return data_CN, data_NC

def process_split(combination, start, end, adata_parts):
    data_CN = []
    data_NC = []

    print('Looking at known ligand-receptor pairs')
    for j in range(start, end):
        print(f"gene pair {combination[j][0]} {combination[j][1]}")
        #print(f"pathway {combination[j][2]}")
        pathway_list.append(combination[j][2])

        # Iterate over all regions
        for region_key, adata_region in adata_parts.items():
            print(f"Processing {region_key}")

            results_CN = lig_rec_plot_normal_cancer_CN(
                adata_region,
                cell_subtype=cell_idx,
                impute=True,
                lig_rec_pair=combination[j],
                all_genes=all_genes
            )
            results_NC = lig_rec_plot_normal_cancer_NC(
                adata_region,
                cell_subtype=cell_idx,
                impute=True,
                lig_rec_pair=combination[j],
                all_genes=all_genes
            )

            data_CN.append([
                combination[j][0], combination[j][1],
                results_CN['normal_receptor_list'].tolist(), results_CN['cancer_ligand_list'].tolist(),
                results_CN['normal_receptor_list_uc'], results_CN['cancer_ligand_list_uc'],
                #results_CN['mean_cancer_ligand'],results_CN['mean_normal_receptor']
                results_CN['n_con_cancer'], results_CN['n_con_normal']
            ])
            data_NC.append([
                combination[j][0], combination[j][1],
                results_NC['normal_ligand_list'].tolist(), results_NC['cancer_receptor_list'].tolist(),
                results_NC['normal_ligand_list_uc'], results_NC['cancer_receptor_list_uc'],
                #results_NC['mean_normal_ligand'],results_NC['mean_cancer_receptor'],
                results_NC['n_con_cancer'], results_NC['n_con_normal']
            ])
    return data_CN, data_NC

def concatenate_lists(group, columns):
    concatenated = {col: [] for col in columns if col not in ['n_con_cancer', 'n_con_normal']}
    sums = {col: 0 for col in ['n_con_cancer', 'n_con_normal']}
    for col in columns:
        if col in ['n_con_cancer', 'n_con_normal']:
            sums[col] = group[col].sum()
        else:
            for item in group[col]:
                if isinstance(item, list):
                    concatenated[col].extend(item)
                else:
                    concatenated[col].append(item)
    concatenated.update(sums)
    return concatenated

def combine_regions(df, columns_to_combine):
    combined_data = []
    
    # Ensure the DataFrame is sorted by 'Ligand', 'Receptor', 'Pathway'
    df_sorted = df.sort_values(by=['Ligand', 'Receptor', 'Pathway']).reset_index(drop=True)
    
    # Group by 'Ligand', 'Receptor', 'Pathway' and aggregate every 9 rows
    num_regions = 9
    num_pairs = len(df) // num_regions
    
    for i in range(num_pairs):
        group = df_sorted.iloc[i*num_regions : (i+1)*num_regions]
        concatenated = concatenate_lists(group, columns_to_combine)
        combined_row = [group.iloc[0]['Ligand'], group.iloc[0]['Receptor'], group.iloc[0]['Pathway']] + [concatenated[col] for col in columns_to_combine]
        combined_data.append(combined_row)
    
    combined_columns = ['Ligand', 'Receptor', 'Pathway'] + columns_to_combine
    return pd.DataFrame(combined_data, columns=combined_columns)

# Define the columns to be concatenated
columns_to_combine_CN = [
    "normal_receptor_list", "cancer_ligand_list",
    "normal_receptor_list_uc", "cancer_ligand_list_uc",
    "n_con_cancer", "n_con_normal"
]

columns_to_combine_NC = [
    "cancer_receptor_list", "normal_ligand_list",
    "cancer_receptor_list_uc", "normal_ligand_list_uc",
    "n_con_cancer", "n_con_normal"
]

if (ligrec_known):
     data_CN1, data_NC1 =process_split_known_ligrec(combination, start_index, end_index, adata_parts)
     df_CN = pd.DataFrame(data_CN1, columns=["Ligand", "Receptor", "Pathway","normal_receptor_list","cancer_ligand_list", "normal_receptor_list_uc","cancer_ligand_list_uc","n_con_cancer","n_con_normal"])
     df_NC = pd.DataFrame(data_NC1, columns=["Ligand", "Receptor", "Pathway","normal_ligand_list","cancer_receptor_list", "normal_ligand_list_uc","cancer_receptor_list_uc","n_con_cancer","n_con_normal"])
     combined_df_CN = combine_regions(df_CN, columns_to_combine_CN)
     combined_df_NC = combine_regions(df_NC, columns_to_combine_NC)
if (all_genes):
     data_CN1, data_NC1 =process_split(combination, start_index, end_index, adata_parts)
     df_CN = pd.DataFrame(data_CN1, columns=["Ligand", "Receptor", "normal_receptor_list","cancer_ligand_list", "normal_receptor_list_uc","cancer_ligand_list_uc","n_con_cancer","n_con_normal"])
     df_NC = pd.DataFrame(data_NC1, columns=["Ligand", "Receptor", "normal_ligand_list","cancer_receptor_list", "normal_ligand_list_uc","cancer_receptor_list_uc","n_con_cancer","n_con_normal"])

print('df_CN.shape')
print(df_CN.shape)

print('df_NC.shape')
print(df_NC.shape)

print('combined_df_CN.shape')
print(combined_df_CN.shape)

print('combined_df_NC.shape')
print(combined_df_NC.shape)

print('type(combined_df_CN[cancer_ligand_list])')
print(type(combined_df_CN['cancer_ligand_list']))

print('combined_df_CN[cancer_ligand_list]')
print(combined_df_CN['cancer_ligand_list'])

print('len(combined_df_CN[cancer_ligand_list])')
print(len(combined_df_CN['cancer_ligand_list']))

print('combined_df_CN[n_con_cancer]')
print(combined_df_CN['n_con_cancer'])

print('combined_df_CN[n_con_normal]')
print(combined_df_CN['n_con_normal'])

print('combined_df_CN[cancer_ligand_list][0])')
print(combined_df_CN['cancer_ligand_list'][0])

print('len(combined_df_CN[cancer_ligand_list][0])')
print(len(combined_df_CN['cancer_ligand_list'][0]))

print('len(combined_df_CN[normal_receptor_list][0])')
print(len(combined_df_CN['normal_receptor_list'][0]))

print('Note: Not applying eval')

print('combined_df_NC.n_con_cancer:\n')
print(combined_df_NC.n_con_cancer)

print('combined_df_CN.n_con_normal:\n')
print(combined_df_CN.n_con_normal)

print('len(df_NC[normal_ligand_list])')
print(len(combined_df_NC['normal_ligand_list']))

print('len(df_NC[cancer_receptor_list])')
print(len(combined_df_NC['cancer_receptor_list']))

#print(len(Endothelial_CN['normal_receptor_list'][0]))

print('type(combined_df_NC[cancer_receptor_list][0]):\n')
print(type(combined_df_NC['cancer_receptor_list'][0]))

corr_coef=correlation_function

print('corr_coef:\n')
print(corr_coef)

corr_connected_list_CN=[]
sig_corr_connected_list_CN=[]
mean_cancer_ligand_CN=[]
mean_normal_receptor_CN=[]

for i in range(len(combined_df_CN)):
    corr_connected_list_CN.append(round(corr_coef(combined_df_CN['cancer_ligand_list'][i], combined_df_CN['normal_receptor_list'][i])[0],4))
    sig_corr_connected_list_CN.append(round(corr_coef(combined_df_CN['cancer_ligand_list'][i], combined_df_CN['normal_receptor_list'][i])[1],4))
    mean_cancer_ligand_CN.append(np.mean(combined_df_CN['cancer_ligand_list'][i]))
    mean_normal_receptor_CN.append(np.mean(combined_df_CN['normal_receptor_list'][i]))

combined_df_CN['corr_connected']=corr_connected_list_CN
combined_df_CN['sig_corr_connected']=sig_corr_connected_list_CN
combined_df_CN['mean_cancer_ligand']=mean_cancer_ligand_CN
combined_df_CN['mean_normal_receptor']=mean_normal_receptor_CN

corr_connected_list_NC=[]
sig_corr_connected_list_NC=[]
mean_normal_ligand_NC=[]
mean_cancer_receptor_NC=[]

for i in range(len(combined_df_NC)):
    corr_connected_list_NC.append(round(corr_coef(combined_df_NC['normal_ligand_list'][i], combined_df_NC['cancer_receptor_list'][i])[0],4))
    sig_corr_connected_list_NC.append(round(corr_coef(combined_df_NC['normal_ligand_list'][i], combined_df_NC['cancer_receptor_list'][i])[1],4))
    mean_normal_ligand_NC.append(np.mean(combined_df_NC['normal_ligand_list'][i]))
    mean_cancer_receptor_NC.append(np.mean(combined_df_NC['cancer_receptor_list'][i]))

combined_df_NC['corr_connected']=corr_connected_list_NC
combined_df_NC['sig_corr_connected']=sig_corr_connected_list_NC
combined_df_NC['mean_normal_ligand']=mean_normal_ligand_NC
combined_df_NC['mean_cancer_receptor']=mean_cancer_receptor_NC
        
subset_label = args.subset_type if args.subset else "full"

combined_df_CN.to_csv(
    path + Normal_cells[cell_idx] + '_imgs/' +
    Normal_cells[cell_idx] + '_' + ct_level + '_' + version +
    '_CN__main_' + all_flag + '_part_' + str(split_part) +
    '_subset_' + subset_label + '_complete_connected.csv',
    index=False
)

combined_df_NC.to_csv(
    path + Normal_cells[cell_idx] + '_imgs/' +
    Normal_cells[cell_idx] + '_' + ct_level + '_' + version +
    '_NC__main_' + all_flag + '_part_' + str(split_part) +
    '_subset_' + subset_label + '_complete_connected.csv',
    index=False
)

#if subset:
#    combined_df_CN.to_csv(path+Normal_cells[cell_idx]+'_imgs/'+Normal_cells[cell_idx]+'_'+ct_level+ '_'+ version+'_CN_'+'_main'+'_'+all_flag+'_'+'part_'+str(split_part)+'_subset_'+subset_type+'_complete_connected.csv', index=False)
#    combined_df_NC.to_csv(path+Normal_cells[cell_idx]+'_imgs/'+Normal_cells[cell_idx]+'_'+ct_level+ '_'+ version+'_NC_'+'_main'+'_'+all_flag+'_'+'part_'+str(split_part)+'_subset_'+subset_type+'_complete_connected.csv', index=False)
#else:
#    combined_df_CN.to_csv(path+Normal_cells[cell_idx]+'_imgs/'+Normal_cells[cell_idx]+'_'+ct_level+ '_'+ version+'_CN_'+'_main'+'_'+all_flag+'_'+'part_'+str(split_part)+'_subset_'+subset+'_complete_connected.csv', index=False)
#    combined_df_NC.to_csv(path+Normal_cells[cell_idx]+'_imgs/'+Normal_cells[cell_idx]+'_'+ct_level+ '_'+ version+'_NC_'+'_main'+'_'+all_flag+'_'+'part_'+str(split_part)+'_subset_'+subset+'_complete_connected.csv', index=False)

print('Finished run successfully')

# Record end time
end_time = time.time()

# Calculate elapsed time
elapsed_time = end_time - start_time
print(f"Elapsed time: {elapsed_time} seconds")


