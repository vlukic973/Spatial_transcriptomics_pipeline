print('Importing packages')

import argparse
import pandas as pd
import os
import sys
import time
import anndata as ad
import numpy as np
import glob

print('Imported packages successfully')

start_time=time.time()

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Aggregate data files and remap cell annotation')
parser.add_argument('--path', required=True, help='Path to the sc and st data')
parser.add_argument('--prefix', required=True, help='Prefix to the cell and gene metadata for spatial')
parser.add_argument('--version', required=True, help='Prefix to the cell and gene metadata for spatial')

args = parser.parse_args()

path = args.path
prefix = args.prefix
version = args.version

# Use glob to get all filenames matching the pattern
#file_pattern = os.path.join("_split.h5ad")

os.chdir(path)
pattern = '*_split.h5ad'

files = sorted(glob.glob(pattern))

print('files:\n')
print(files)

# Initialize a list to hold the AnnData objects
adata_list = []

# Loop through the files and read each AnnData object
for file in files:
    if os.path.getsize(file) > 0:  # Ensure the file is not empty
        print(f"Reading {file}")
        adata = ad.read_h5ad(file)
        if adata.shape[0] > 0:  # Ensure the AnnData object is not empty
            adata_list.append(adata)
        else:
            print(f"Warning: {file} contains no cells. Skipping.")
    else:
        print(f"Warning: {file} is empty. Skipping.")

# Check if we have any valid AnnData objects to concatenate
if len(adata_list) == 0:
    raise ValueError("No valid AnnData objects found to concatenate.")

# Concatenate all the AnnData objects
adata_all = ad.concat(adata_list, axis=0, join='outer')

print('adata_all:\n')
print(adata_all)

print(np.unique(adata_all.obs['tangram_ct']))

cell_type_mapping = {
    'Endothelial': 'Endothelial',
    'Fibroblast': 'Fibroblast',
    'Hepatocyte': 'Epithelial',  # Hepatocytes as epithelial cells
    'Malignant': 'Epithelial',  # Malignant cells assumed to be epithelial
    'Myeloid': 'Myeloid',
    'Pericyte': 'Unknown',  # Or keep as 'Pericyte' if needed separately
    'Plasma': 'Lymphocyte',
    'T': 'Lymphocyte',
    'Unknown': 'Unknown',  # Excluding unknowns
    'myofibroblast': 'Fibroblast',
    'B': 'Lymphocyte',
    'Epithelial': 'Epithelial',
    'NK': 'Lymphocyte',
    'cycling': 'Unknown' 
}

# Function to map cell types
def map_cell_types(cell_type):
    return cell_type_mapping.get(cell_type, 'Unknown')

#example_cell_types = ['Hepatocyte', 'Malignant', 'Endothelial', 'Unknown', 'Pericyte', 'cycling']

cell_types_to_remap=adata_all.obs['tangram_ct']

mapped_cell_types = [map_cell_types(ct) for ct in cell_types_to_remap]

print('np.unique(mapped_cell_types)')
print(np.unique(mapped_cell_types))  # Output based on mapping rules

adata_all.obs['ct_level0']=mapped_cell_types

# Filter out 'Unknown' cell types
adata_filtered = adata_all[adata_all.obs['ct_level0'] != 'Unknown']

# Save the AnnData object

output_file = os.path.join(path,"adata_"+version+".h5ad")
adata_filtered.write_h5ad(output_file)
print(f"Aggregated AnnData object saved to {output_file}")




