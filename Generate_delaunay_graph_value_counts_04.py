print('Running GR_sample_generate_delaunay_graph script')

import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import anndata as ad
from pathlib import Path
import scipy
from matplotlib.pyplot import rc_context
import argparse
import squidpy as sq
sc.logging.print_header()

# Create the argument parser
parser = argparse.ArgumentParser(description='Generate delaunay graph')

# Add the --path argument
parser.add_argument('--path', type=str, help='Path to files')
parser.add_argument('--ct_level', type=str, help='Cell-type annotation level')
parser.add_argument('--version', type=str, help='Run version')
parser.add_argument('--prefix', type=str, help='Prefix')

# Parse the command-line arguments
args = parser.parse_args()

path = args.path
ct_level = args.ct_level
version = args.version
prefix = args.prefix

# Use the value of the --path argument in your script
# For example:
print(f'The specified path is: {path}')

file='SCVI_imputed_anndata_GR_'+version+'.h5ad'

print('Reading in imputed anndata matrix')
adata=ad.read_h5ad(path+file)

print('adata:\n')
print(adata)

# Stack the arrays along a new axis
xy_coords = np.stack((adata.obs['center_x'],adata.obs['center_y']), axis=1)

print('xy_coords:\n')
print(xy_coords)

adata.obsm['spatial']=xy_coords

#print('adata.obs.index')
#print(adata.obs.index)

print(adata.obsm['spatial'])

print('Doing sq.gr.spatial_neighbors')

sq.gr.spatial_neighbors(adata, coord_type="generic", radius=(0, 40), delaunay=True)
adata.uns['spatial_neighbors']['params']['radius']=list(adata.uns['spatial_neighbors']['params']['radius'])

print('Writing h5ad object with Delaunay graph')
print(path+'anndata_GR_with_delaunay_40_microns_'+ct_level+'_'+version+'.h5ad')
adata.write_h5ad(path+'anndata_GR_with_delaunay_40_microns_'+ct_level+'_'+version+'.h5ad')

print("Here adata.obs.columns")
print(adata.obs.columns)

#if hasattr(adata.obs, 'ct_level0'):
    # Create the scatter plot
fig, ax = plt.subplots(figsize=(15, 10))
sq.pl.spatial_scatter(
adata,
color=ct_level,
connectivity_key="spatial_connectivities",
edges_color="black",
shape=None,
edges_width=2,
size=10,
ax=ax
)

# Set the background color to white
fig.set_facecolor('white')

# Save the figure
fig.savefig(path+'Delaunay_graph_connected_'+ct_level+'_'+version+'.png')

plt.rcParams['axes.facecolor']='white'
plt.rcParams['savefig.facecolor']='white'
plt.figure(figsize=(8, 8))
adata.obs[ct_level].value_counts().plot(kind='bar')
plt.subplots_adjust(bottom=0.3)
plt.savefig(path+'Value_counts_histogram_'+ct_level+'_'+version+'.png')
