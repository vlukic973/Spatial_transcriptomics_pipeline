This document explains how to use the pipeline to get candidate direct-contact ligand-receptor pairs from an ovarian cancer cell sample.

1. Download publicly available cancer samples from Vizgen website, as well as a single cell reference for a particular cancer type.

The publicly available cancer samples will look like e.g. 

- HumanOvarianCancerPatient1_cell_by_gene.csv
- HumanOvarianCancerPatient1_cell_metadata.csv

for a human ovarian cancer sample.

The single-cell reference will look like e.g.:

OCA_UNB_10X_E-MTAB-8107.h5ad

Make the files into an anndata (adata) object.

2. Run Get_celltype_annotation_tangram_01.py which splits up the anndata object, does the cell type annotation from using a single cell reference and mapping it over (necessary because tangram runs into memory issues if we are putting in the entire adata object if it has more than around 100,000 cells).

3. Run Aggregate_adata_remap_cells_02.py. This aggregates the annotated split files and does a remapping into the coarsest annotations in the categories of Epithelial, Myeloid, Endothelial, Lymphocyte
and Fibroblast.

4. Run Pre_process_data_cluster_03.py. This script imputes the zero-valued expression data using SCVI. 

5. Run Generate_delaunay_graph_value_counts_04.py. This makes a Delaunay graph for the anndata object, for which we use a threshold of 40 microns to say that cells are in direct contact.

6. Run Submit_ligand_receptor_celltype_scripts_05.py. This script takes as input either 1) all pairwise gene combinations, 2) Only DE genes or 3) Known ligand-receptor pairs (obtained from the CellChat database). We have mainly focused on the known-ligand receptor pairs in this project. The script then activates another script based on whether 1), 2) or 3) are chosen, which then calls 


Perform the ligand-receptor analysis (using Ligrec_split_run_gene_pair_stats_ctlevel_connected_unconnected_all_dirs_expression_public_v3.py), subsetting the list to known direct-
contact ligand-receptor pairs from CellChat. Calculate the correlation between each ligand-receptor pair, for cancer cells in contact with normal cell subtypes and vice versa. This involves splitting up the original anndata object into 9 parts (necessary for the liver and lung, as number of cells is very large)

7. Run Known_all_ligrec_pairs_stats_cluster_v3.py to perform a significance test (using z-test)
with additional requirement that the expression of a particular ligand-receptor pair has to be greater than the mean.

8. Assess quality of ligand-receptor pairs output by the method, using a Machine-learning method (specifically, random forest). The expression of the connected cells between cancer and normal subtypes is used, for each particular ligand-receptor pair. We use a label of 1 if a particular ligand-receptor pair is predicted by a particular method (currently it is our method versus Cellchat).

The script was run using the following command:

python ML_metrics_CellChat_Home_LIANA_subset.py --cell_subtype "Endothelial" --ct_level "ct_level0" --sample "Ovarian" --version "v3" --path "/path/to/file/HumanOvarianCancerPatient1/"

We noted later that this method of assessing ligand-receptor pairs was biased, as it relied on the gene expression data from connected cells only, and this is the data used by our method.
