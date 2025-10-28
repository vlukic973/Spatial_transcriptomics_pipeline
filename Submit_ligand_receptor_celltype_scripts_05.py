print("Importing packages including anndata")

import math
import argparse
import numpy as np
import subprocess
import pandas as pd
import anndata as ad
import os

def divisorGenerator(n):
    divisors = []
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divisors.append(i)
            if i != n // i:  # Only add non-square divisors
                divisors.append(n // i)
    divisors.sort()  # Ensure divisors are in ascending order
    return divisors

def findCorrectDivisors(number):
    divisors = divisorGenerator(number)
    total_divisors = len(divisors)
    mid = total_divisors // 2
    
    # Check the middle two divisors first
    div1, div2 = divisors[mid - 1], divisors[mid]
    if div1 * div2 == number:
        return div1, div2
    
    # If the middle two divisors didn't work, check other pairs of divisors
    for i in range(mid - 1, -1, -1):
        div1 = divisors[i]
        div2 = number // div1
        if div1 * div2 == number:
            return div1, div2

    return None  # If no suitable divisors are found

def count_lines_in_file(filename):
    with open(filename, 'r') as file:
        return sum(1 for line in file)

def main():
    parser = argparse.ArgumentParser(description="Generate divisors for a given number or file with line count.")
    parser.add_argument("--path", type=str, help="The base path for the file.")
    parser.add_argument("--cell_subtype", type=str, help="The cell subtype.")
    parser.add_argument("--input", type=str, help="The number or text file with line count.")
    parser.add_argument("--ct_level", type=str, help="Cell type annotation level")
    parser.add_argument("--version", type=str, help="Version number")
    parser.add_argument("--corr_coef", type=str, help="Pearson or Spearman")
    parser.add_argument("--sample", type=str, help="Which sample")
    parser.add_argument("--all_genes", action='store_true', help='Consider all pairwise genes')
    parser.add_argument("--ligrec_known", action='store_true', help='Consider known ligrec pairs')
    parser.add_argument("--DE_genes", action='store_true', help='Consider DE genes')
    parser.add_argument("--subset", action='store_true', help='Consider a subset of the original data')
    parser.add_argument("--subset_type", type=str, choices=['60p', '70p', '80p'], help='Which subset to use: 60p, 70p, or 80p')

    args = parser.parse_args()

    all_genes = False
    ligrec_known = False
    DE_genes = False
    subset = False

    if args.all_genes:
        print("Considering all genes")
        all_genes = True
    if args.ligrec_known:
        print("Considering all known ligrec genes")
        ligrec_known = True
    if args.DE_genes:
        print("Considering DE-genes")
        DE_genes = True
    if args.subset:
        print("Considering data subset")
        subset = True

    if args.input:
        if args.input.isdigit():
            number = int(args.input)
        else:
            number = count_lines_in_file(args.input)
    else:
        if all_genes:
            print('Reading in all gene pair combinations')
            filename = "/path/to/file/all_gene_pairs_TNBC.txt"
            number = count_lines_in_file(filename)
        if ligrec_known:
            print("Considering known ligrec_pairs")
            #filename = "{}/gene_pairs_known_ligrec.txt".format(args.path)
            #print('filename:\n')
            #print(filename)
            print("Getting known ligand-receptor pairs")
            adata = ad.read_h5ad(os.path.join(args.path,'anndata_GR_with_delaunay_40_microns_'+args.ct_level+'_'+args.version+'.h5ad'))

            ligand_receptor_df=pd.read_csv('/path/to/file/'+'ligand_receptor_pairs.txt',sep=r'\s+',on_bad_lines='skip')

            # Split ligands and receptors if they contain entries separated by '_'
            ligand_receptor_df[['ligand1', 'ligand2']] = ligand_receptor_df['ligand'].str.split('_', n=1, expand=True)
            ligand_receptor_df[['receptor1', 'receptor2']] = ligand_receptor_df['receptor'].str.split('_', n=1, expand=True)

            ligand_receptor_df_subset=ligand_receptor_df[['ligand1','ligand2','receptor1','receptor2']]

            ligand_receptor_df_subset = ligand_receptor_df_subset.applymap(lambda x: x.upper() if isinstance(x, str) else x)

            panel_genes = pd.Series(adata.var_names)

            panel_genes = panel_genes.str.upper()

            # Initialize an empty list to store valid pairs
            valid_pairs = []

            # Iterate through rows
            for index, row in ligand_receptor_df_subset.iterrows():
                # Iterate through combinations
                for ligand, receptor in [(row['ligand1'], row['receptor1']),
                             (row['ligand1'], row['receptor2']),
                             (row['ligand2'], row['receptor1']),
                             (row['ligand2'], row['receptor2'])]:
                # Check if both ligand and receptor are in panel_genes
                    if ligand in list(panel_genes) and receptor in list(panel_genes):
                        valid_pairs.append((ligand, receptor))
            #pathway_list.append(pathway)

            print('valid_pairs:\n')
            print(valid_pairs)

            with open(args.path+"gene_pairs_cellchatdb_direct_contact.txt", "w") as file:
                for pair in valid_pairs:
                    file.write(" ".join(pair) + "\n")

            filename = "{}/gene_pairs_cellchatdb_direct_contact.txt".format(args.path)
            number = count_lines_in_file(filename)
        if DE_genes:
            print('Reading in DE gene pair combinations')
            filename = "{}/{}_gene_pairs_CN_{}_{}_{}_cat.txt".format(args.path, args.cell_subtype, args.ct_level, args.version, args.sample)
            number = count_lines_in_file(filename)

    print("Input number: {}".format(number))
    divisors = findCorrectDivisors(number)
    if divisors:
        print(f"The correct divisors for {number} are: {divisors[0]}, {divisors[1]}")
    else:
        print("No suitable divisors found. Exiting...")
        exit()
    div1, div2 = divisors

#    if ligrec_known:

    if all_genes:
        print('Executing sbatch for all genes using run_script_DE_parallel_ctlevel_all_vs_all.sh')
        command = [
            "sbatch",
            "-e", "my_job_error_test_all_finer_full.log",
            "--array=1-{}".format(div1),
            "run_script_DE_parallel_ctlevel_all_vs_all.sh",
            str(div1), str(div2), args.path, args.cell_subtype,
            args.ct_level, args.version, args.corr_coef, args.sample
        ]
    elif ligrec_known:
        print('Executing sbatch for all genes using run_script_DE_parallel_ctlevel_known_ligrec.sh')
        command = [
            "sbatch",
            "-e", "my_job_error_test_all_coarse_full.log",
            "--array=1-{}".format(div1),
            "run_script_ctlevel_known_ligrec.sh",
            str(div1), str(div2), args.path, args.cell_subtype,
            args.ct_level, args.version, args.corr_coef, args.sample, args.subset_type
        ]
    elif DE_genes:
        print('Executing sbatch for DE genes')
        command = [
            "sbatch",
            "-e", "my_job_error_test_" + args.cell_subtype + "_finer_full.log",
            "--array=1-{}".format(div1),
            "run_script_DE_parallel_ctlevel_all.sh",
            str(div1), str(div2), args.path, args.cell_subtype,
            args.ct_level, args.version, args.corr_coef, args.sample
        ]

    print('Check squeue')
    # Execute the command
    print('Executing: \n')
    print(command)
    subprocess.run(command)

if __name__ == "__main__":
    main()
