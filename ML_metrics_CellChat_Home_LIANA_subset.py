import numpy as np
import pandas as pd
import ast
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, f1_score
import argparse
import os
import matplotlib.pyplot as plt
import csv

def write_results(patient, cell_type, subset_type, f1_home, f1_cellchat, f1_random_tuple, f1_liana_spec, f1_liana_mag, output_file="ML_metrics_results.csv"):
    header = ["Patient", "cell_type", "subset_type","F1_home", "F1_CellChat", "F1_random", "F1_LIANA_spec", "F1_LIANA_mag"]
    file_exists = os.path.isfile(output_file)
    
    with open(output_file, mode='a', newline='') as file:
        writer = csv.writer(file)
        
        # Write header only if the file does not exist (first run)
        if not file_exists:
            writer.writerow(header)
        
        # Append the new row with F1_random as a tuple string
        writer.writerow([patient, cell_type, subset_type, f1_home if f1_home is not None else "NA", f1_cellchat, f1_random_tuple, f1_liana_spec, f1_liana_mag])


parser = argparse.ArgumentParser(description='Perform ML validation of Cellchat and Home method ligand-receptor pairs')

parser.add_argument('--cell_subtype', type=str, help='Normal cell subtype')
parser.add_argument('--ct_level', type=str, required=True, help='Cell type annotation')
parser.add_argument('--sample',type=str,help="Which sample")
parser.add_argument('--version', type=str, required=True, help='Run version')
parser.add_argument('--path', type=str, required=True, help='path to home data')
parser.add_argument('--subset_type', type=str, choices=['60p', '70p', '80p'], help='Which subset to use: 60p, 70p, or 80p')

args = parser.parse_args()

cell_subtype = args.cell_subtype
ct_level=args.ct_level
sample=args.sample
version=args.version
path=args.path
subset_type=args.subset_type

# Handle subset_type only if provided
subset_type = args.subset_type if args.subset_type is not None else None
#subset_type = args.subset_type if args.subset_type else "Complete"

# Example usage
if subset_type:
    print(f"Using subset type: {subset_type}")
else:
    print("No subset type provided.")
    #subset_type

### Home method NC

# Function to save significant pairs to CSV
def save_significant_pairs(pairs, method, cell_subtype, sample, version, path):
    # Create the directory if it does not exist
    output_dir = os.path.join(path, "results", sample, cell_subtype, version)
    os.makedirs(output_dir, exist_ok=True)

    # Define the output file path
    output_file = os.path.join(output_dir, f"significant_pairs_{method}.csv")

    # Save the pairs to a CSV file
    pairs_df = pd.DataFrame(pairs, columns=['Ligand', 'Receptor'])
    pairs_df.to_csv(output_file, index=False)
    print(f"Saved {method} significant pairs to {output_file}")

# Extract significant ligand-receptor pairs
def extract_significant_pairs(df):
    significant_pairs = df[df['significant'] == 1][['Ligand', 'Receptor']]
    return [tuple(x) for x in significant_pairs.to_records(index=False)]

# Jaccard index function
def jaccard_index(set1, set2):
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    return intersection / union if union > 0 else 0

# Replace 'nan' strings with 'np.nan' before applying eval
def is_nan_list(lst):
    # Check if the list is indeed a list and all its elements are NaN
    return isinstance(lst, list) and all(np.isnan(x) for x in lst)

#print('Reading in connected_full_NC')

connected_full_NC=pd.read_csv(path+cell_subtype+'_'+sample+'_concatenated_'+'NC_'+ct_level+'_'+version+'.csv')

# Drop rows where 'corr_connected' column has NaN values
connected_full_NC = connected_full_NC.dropna(subset=['corr_connected'])

# Endothelial_ovarian_bh_correction_pval_mean_filtered_CN_ct_level0_0.05_v3_known_ligrec.csv

#print('Reading in sig_NC')

#file_path = path + cell_subtype + '_' + sample + '_' + 'bh_correction_pval_mean_filtered_NC_' + ct_level + '_0.05_' + version + '_known_ligrec.csv'
if subset_type:
    file_path = path + cell_subtype + '_' + sample + '_' + 'bh_correction_pval_mean_filtered_NC_' + ct_level+'_subset'+subset_type+ '_0.05_' + version + '_known_ligrec.csv'
    print('Read in NC home predictions from '+ subset_type+ ' dataset')
else:
    file_path = path + cell_subtype + '_' + sample + '_' + 'bh_correction_pval_mean_filtered_NC_' + ct_level + '_0.05_' + version + '_known_ligrec.csv'
    print('Read in NC home predictions from full dataset')

if os.path.exists(file_path):
    sig_NC = pd.read_csv(file_path)
else:
    print(f"File does not exist: {file_path}")
    sig_NC = pd.DataFrame()

#sig_NC=pd.read_csv(path+cell_subtype+'_'+sample+'_'+'bh_correction_pval_mean_filtered_NC_'+ct_level+'_0.05_'+version+'_known_ligrec.csv')

if not sig_NC.empty:
    print("The DataFrame sig_NC is not empty. Continuing analysis")
    sig_NC = sig_NC.drop('Unnamed: 0', axis=1)
    sig_NC = sig_NC.drop_duplicates()
    # Merge df1 with df2 to find common 'Ligand' and 'Receptor' pairs
    merged_NC = connected_full_NC.merge(sig_NC, on=['Ligand', 'Receptor'], how='left', indicator=True)
    # Assign 'significant' value based on the merge result
    connected_full_NC['significant'] = merged_NC['_merge'].apply(lambda x: 1 if x == 'both' else 0)
    connected_full_NC=connected_full_NC.fillna(0)
    #print('Merged home_NC:')
    #print(connected_full_NC[['Ligand', 'Receptor','significant']])
    connected_full_NC_Jaccard_home=connected_full_NC[['Ligand', 'Receptor','significant','corr_connected', 'sig_corr_connected']]
    #print('significant only connected_full_NC_Jaccard_home')
    #print(connected_full_NC_Jaccard_home[connected_full_NC_Jaccard_home['significant']==1].drop_duplicates())
    connected_full_NC=connected_full_NC[['normal_ligand_list','cancer_receptor_list','significant']]
    #print('Here connected_full_NC')
    #print(connected_full_NC)
    connected_full_NC['normal_ligand_list'] = connected_full_NC['normal_ligand_list'].apply(eval)
    connected_full_NC['cancer_receptor_list'] = connected_full_NC['cancer_receptor_list'].apply(eval)
    # Expand the normal_ligand_list into separate columns
    normal_ligand_expanded = connected_full_NC['normal_ligand_list'].apply(lambda x: pd.Series(x))
    normal_ligand_expanded.columns = [f'normal_ligand_{i}' for i in normal_ligand_expanded.columns]
    # Expand the cancer_receptor_list into separate columns
    cancer_receptor_expanded = connected_full_NC['cancer_receptor_list'].apply(lambda x: pd.Series(x))
    cancer_receptor_expanded.columns = [f'cancer_receptor_{i}' for i in cancer_receptor_expanded.columns]
    # Concatenate the expanded columns back with the original DataFrame (dropping the original list columns)
    df_expanded_NC = pd.concat([connected_full_NC.drop(['normal_ligand_list', 'cancer_receptor_list'], axis=1), normal_ligand_expanded, cancer_receptor_expanded], axis=1)
    df_expanded_NC_cols = pd.DataFrame(columns=[f'column_{i}' for i in range(df_expanded_NC.shape[1])])
    # Generate new column names
    new_column_names = [str(i) for i in range(1, len(df_expanded_NC_cols.columns) + 1)]
    df_expanded_NC.columns = new_column_names
    df_expanded_NC_expression=df_expanded_NC.copy()
    df_expanded_NC_expression=pd.concat([connected_full_NC_Jaccard_home[['Ligand', 'Receptor','corr_connected', 'sig_corr_connected']],df_expanded_NC], axis=1)
    #print('df_expanded_NC_expression')
    #print(df_expanded_NC_expression)
    rows_of_interest = df_expanded_NC_expression[df_expanded_NC_expression['1'] == 1]
    # Get the midpoint column
    #print('df_expanded_NC_expression.shape[1]')
    #print(df_expanded_NC_expression.shape[1])
    midpoint = (df_expanded_NC_expression.shape[1]-5) // 2
    midpoint=midpoint+3
    print(midpoint)
    # Step 2: Iterate through these rows and plot ligand vs receptor expression
    for index, row in rows_of_interest.iterrows():
        ligand_expression = row[2:midpoint]  # Columns 2 to 2084/2
        receptor_expression = row[midpoint+1:df_expanded_NC_expression.shape[1]]  # Columns 2084/2 to 2084
        #print('len(ligand_expression)')
        #print(len(ligand_expression))
        #print('len(receptor_expression)')
        #print(len(receptor_expression))
        # Retrieve the correlation and significance values
        corr_connected = row['corr_connected']
        sig_corr_connected = row['sig_corr_connected']
    # Step 3: Plot ligand vs receptor
        plt.figure(figsize=(10, 6))
        plt.scatter(ligand_expression, receptor_expression, label=f"{row['Ligand']} - {row['Receptor']}")
        plt.xlabel(f"Ligand Expression ({row['Ligand']})")
        plt.ylabel(f"Receptor Expression ({row['Receptor']})")
        plt.title(f"Ligand ({row['Ligand']}) vs Receptor ({row['Receptor']}) Expression")
        plt.title(f"Ligand ({row['Ligand']}) vs Receptor ({row['Receptor']}) Expression\n"
              f"Correlation: {corr_connected:.2f}, Significance: {sig_corr_connected:.2e}")
        plt.legend()
        plt.savefig(path+'Home_NC_'+cell_subtype+'_'+row['Ligand']+'_'+row['Receptor']+'_corr_plot.png')
        plt.close()

    # Assign new column names to the dataframedf_expanded_NC.columns = new_column_names
    # Proceed with your operations
else:
    print("The DataFrame sig_NC is empty.")
    df_expanded_NC=pd.DataFrame()

### Home method CN

#connected_full_CN=pd.read_csv(path+cell_subtype+'_'+ct_level+'_'+version+'_CN__main_known_ligrec_genes_part_1_complete_connected.csv')
connected_full_CN=pd.read_csv(path+cell_subtype+'_'+sample+'_concatenated_'+'CN_'+ct_level+'_'+version+'.csv')
print("Successfully read in connected_full_CN")

# Drop rows where 'corr_connected' column has NaN values
connected_full_CN = connected_full_CN.dropna(subset=['corr_connected'])

if subset_type:
#file_path = path + cell_subtype + '_' + sample + '_' + 'bh_correction_pval_mean_filtered_CN_' + ct_level + '_0.05_' + version + '_known_ligrec.csv'
    file_path = path + cell_subtype + '_' + sample + '_' + 'bh_correction_pval_mean_filtered_CN_' + ct_level +'_subset'+subset_type+ '_0.05_' + version + '_known_ligrec.csv'
    print('Read in CN home predictions from '+ subset_type+ ' dataset')
else:
    file_path = path + cell_subtype + '_' + sample + '_' + 'bh_correction_pval_mean_filtered_CN_' + ct_level + '_0.05_' + version + '_known_ligrec.csv'
    print('Read in CN home predictions from full dataset')

if os.path.exists(file_path):
    sig_CN = pd.read_csv(file_path)
else:
    print(f"File does not exist: {file_path}")
    sig_CN = pd.DataFrame()

if not sig_CN.empty:
    sig_CN = sig_CN.drop('Unnamed: 0', axis=1)
    sig_CN = sig_CN.drop_duplicates()
    merged_CN = connected_full_CN.merge(sig_CN, on=['Ligand', 'Receptor'], how='left', indicator=True)
    connected_full_CN['significant'] = merged_CN['_merge'].apply(lambda x: 1 if x == 'both' else 0)
    connected_full_CN=connected_full_CN.fillna(0)
    connected_full_CN_Jaccard_home=connected_full_CN[['Ligand', 'Receptor','significant','corr_connected', 'sig_corr_connected']]
    connected_full_CN=connected_full_CN[['normal_receptor_list','cancer_ligand_list','significant']]
    connected_full_CN['normal_receptor_list'] = connected_full_CN['normal_receptor_list'].apply(eval)
    connected_full_CN['cancer_ligand_list'] = connected_full_CN['cancer_ligand_list'].apply(eval)
    normal_receptor_expanded = connected_full_CN['normal_receptor_list'].apply(lambda x: pd.Series(x))
    normal_receptor_expanded.columns = [f'normal_receptor_{i}' for i in normal_receptor_expanded.columns]
    # Expand the cancer_ligand_list into separate columns
    cancer_ligand_expanded = connected_full_CN['cancer_ligand_list'].apply(lambda x: pd.Series(x))
    cancer_ligand_expanded.columns = [f'cancer_ligand_{i}' for i in cancer_ligand_expanded.columns]
    # Concatenate the expanded columns back with the original DataFrame (dropping the original list columns)
    df_expanded_CN = pd.concat([connected_full_CN.drop(['normal_receptor_list', 'cancer_ligand_list'], axis=1), normal_receptor_expanded, cancer_ligand_expanded], axis=1)
    df_expanded_CN_cols = pd.DataFrame(columns=[f'column_{i}' for i in range(df_expanded_CN.shape[1])])
    # Generate new column names
    new_column_names = [str(i) for i in range(1, len(df_expanded_CN_cols.columns) + 1)]
    # Assign new column names to the dataframe
    df_expanded_CN.columns = new_column_names
    #print('df_expanded_CN.shape:')
    #print(df_expanded_CN.shape)
    df_expanded_CN_expression=df_expanded_CN.copy()
    df_expanded_CN_expression=pd.concat([connected_full_CN_Jaccard_home[['Ligand', 'Receptor','corr_connected', 'sig_corr_connected']],df_expanded_CN], axis=1)
    #print('df_expanded_CN_expression')
    #print(df_expanded_CN_expression)
    rows_of_interest = df_expanded_CN_expression[df_expanded_CN_expression['1'] == 1]
    # Get the midpoint column
    #print('df_expanded_CN_expression.shape[1]')
    #print(df_expanded_CN_expression.shape[1])
    midpoint = (df_expanded_CN_expression.shape[1]-5) // 2
    midpoint=midpoint+3
    #print(midpoint)
    # Step 2: Iterate through these rows and plot ligand vs receptor expression
    for index, row in rows_of_interest.iterrows():
        ligand_expression = row[2:midpoint]  # Columns 2 to 2084/2
        receptor_expression = row[midpoint+1:df_expanded_CN_expression.shape[1]]  # Columns 2084/2 to 2084
        corr_connected = row['corr_connected']
        sig_corr_connected = row['sig_corr_connected']
    # Step 3: Plot ligand vs receptor
        plt.figure(figsize=(10, 6))
        plt.scatter(ligand_expression, receptor_expression, label=f"{row['Ligand']} - {row['Receptor']}")
        plt.xlabel(f"Ligand Expression ({row['Ligand']})")
        plt.ylabel(f"Receptor Expression ({row['Receptor']})")
        plt.title(f"Ligand ({row['Ligand']}) vs Receptor ({row['Receptor']}) Expression")
        plt.title(f"Ligand ({row['Ligand']}) vs Receptor ({row['Receptor']}) Expression\n"
              f"Correlation: {corr_connected:.2f}, Significance: {sig_corr_connected:.2e}")
        plt.legend()
        plt.savefig(path+'Home_CN_'+cell_subtype+'_'+row['Ligand']+'_'+row['Receptor']+'_corr_plot.png')
        plt.close()
else:
    print("The DataFrame sig_CN is empty.")
    df_expanded_CN=pd.DataFrame()

dfs_to_concat = [df for df in [df_expanded_CN, df_expanded_NC] if not df.empty]

# Concatenate only non-empty DataFrames
df_expanded_all = pd.concat(dfs_to_concat, ignore_index=True) if dfs_to_concat else pd.DataFrame()

if not df_expanded_all.empty:
    features = df_expanded_all.drop('1', axis=1)
    target = df_expanded_all['1']
    #print('Here is target.value_counts() for home method')
    print('target.value_counts()')
    print(target.value_counts())
    # Ensure both classes exist before accessing them
    value_counts = target.value_counts()
    if all(cls in value_counts for cls in [0, 1]) and value_counts[1] > 1 and value_counts[0] > 1:
        print('There are enough instances in the classes')

        print('Splitting into train/valid/test sets')
        X_train, X_test, y_train, y_test = train_test_split(features, target, test_size=0.2, random_state=42, stratify=target)

# Initialize the Random Forest classifier
        print('Initialising the Random Forest classifier')
        rf_model = RandomForestClassifier(n_estimators=100, random_state=42)

# Train the model on the training data
        print('Train the model on the training data')
        rf_model.fit(X_train, y_train)

    # Make predictions on the test set
        y_test_pred = rf_model.predict(X_test)

    # Calculate accuracy
        test_accuracy = accuracy_score(y_test, y_test_pred)
        print(f'Home test Accuracy: {test_accuracy:.4f}')
        # Print classification report
        print("Home classification report:\n", classification_report(y_test, y_test_pred))

        f1_home=f1_score(y_test, y_test_pred, pos_label=1)
        print('f1 score for home is '+ str(f1_home))

    # Print confusion matrix
        #print("Home confusion matrix:\n", confusion_matrix(y_test, y_test_pred))
    else:
        print('There are not enough instances in the classes')
        f1_home=None
else:
    print('df_expanded_all dataframe is empty so could not do ML part using home method')
    f1_home=None

print('Finished checking home part')

#################
####CELLCHAT#####
#################

print('Checking CellChat results')

#NC_Cellchat=pd.read_csv(path+sample+'_CellChat_ligand_receptor_pair_counts_'+cell_subtype+'_to_Epithelial_only_all_present.csv')
#CN_Cellchat=pd.read_csv(path+sample+'_CellChat_ligand_receptor_pair_counts_Epithelial_to_'+cell_subtype+'_only_all_present.csv')

if subset_type:
    NC_Cellchat=pd.read_csv(path+sample+'_CellChat_ligand_receptor_pair_counts_'+cell_subtype+'_to_Epithelial_'+subset_type+'_only_all_present.csv')
    CN_Cellchat=pd.read_csv(path+sample+'_CellChat_ligand_receptor_pair_counts_Epithelial_to_'+cell_subtype+'_'+subset_type+'_only_all_present.csv')
    print('Read in cellchat subset '+subset_type+' predictions')
else:
    NC_Cellchat=pd.read_csv(path+sample+'_CellChat_ligand_receptor_pair_counts_'+cell_subtype+'_to_Epithelial_only_all_present.csv')
    CN_Cellchat=pd.read_csv(path+sample+'_CellChat_ligand_receptor_pair_counts_Epithelial_to_'+cell_subtype+'_only_all_present.csv')
    print('Read in full cellchat predictions')

#print(sample+'_CellChat_ligand_receptor_pair_counts_'+cell_subtype+'_to_Epithelial_only_all_present.csv')
#print(sample+'_CellChat_ligand_receptor_pair_counts_Epithelial_to_'+cell_subtype+'_only_all_present.csv')

connected_full_NC=pd.read_csv(path+cell_subtype+'_'+sample+'_concatenated_'+'NC_'+ct_level+'_'+version+'.csv')

print('Dropping NAs in cellchat results NC')
connected_full_NC = connected_full_NC.dropna(subset=['corr_connected'])

merged_cellchat_NC = connected_full_NC.merge(NC_Cellchat, on=['Ligand', 'Receptor'], how='left', indicator=True)
connected_full_NC['significant'] = merged_cellchat_NC['_merge'].apply(lambda x: 1 if x == 'both' else 0)
#print('Merged cellchat_NC:')
#print(connected_full_NC[['Ligand', 'Receptor','significant']])
connected_full_NC_Jaccard_cellchat=connected_full_NC[['Ligand', 'Receptor','significant']]
#print('Cellchat NC significant only')
#print(connected_full_NC_Jaccard_cellchat[connected_full_NC_Jaccard_cellchat['significant']==1].drop_duplicates())
connected_full_NC=connected_full_NC[['normal_ligand_list','cancer_receptor_list','significant']]
connected_full_NC['normal_ligand_list'] = connected_full_NC['normal_ligand_list'].apply(eval)
connected_full_NC['cancer_receptor_list'] = connected_full_NC['cancer_receptor_list'].apply(eval)
normal_ligand_expanded = connected_full_NC['normal_ligand_list'].apply(lambda x: pd.Series(x))
normal_ligand_expanded.columns = [f'normal_ligand_{i}' for i in normal_ligand_expanded.columns]
cancer_receptor_expanded = connected_full_NC['cancer_receptor_list'].apply(lambda x: pd.Series(x))
cancer_receptor_expanded.columns = [f'cancer_receptor_{i}' for i in cancer_receptor_expanded.columns]
df_expanded_NC = pd.concat([connected_full_NC.drop(['normal_ligand_list', 'cancer_receptor_list'], axis=1), normal_ligand_expanded, cancer_receptor_expanded], axis=1)

df_expanded_NC_cols = pd.DataFrame(columns=[f'column_{i}' for i in range(df_expanded_NC.shape[1])])
# Generate new column names
new_column_names = [str(i) for i in range(1, len(df_expanded_NC_cols.columns) + 1)]
# Assign new column names to the dataframe
df_expanded_NC.columns = new_column_names

#connected_full_CN=pd.read_csv(path+'Lymphocyte_ct_level0_raw_counts_v2_CN__main_known_ligrec_genes_part_1_complete_connected.csv')
#connected_full_CN=pd.read_csv(path+cell_subtype+'_'+ct_level+'_'+version+'_CN__main_known_ligrec_genes_part_1_complete_connected.csv')
connected_full_CN=pd.read_csv(path+cell_subtype+'_'+sample+'_concatenated_'+'CN_'+ct_level+'_'+version+'.csv')
#print('Dropping NAs in cellchat results')
connected_full_CN = connected_full_CN.dropna(subset=['corr_connected'])
merged_cellchat_CN = connected_full_CN.merge(CN_Cellchat, on=['Ligand', 'Receptor'], how='left', indicator=True)
connected_full_CN['significant'] = merged_cellchat_CN['_merge'].apply(lambda x: 1 if x == 'both' else 0)
#print('Merged cellchat_CN:')
#print(connected_full_CN[['Ligand', 'Receptor','significant']])
connected_full_CN_Jaccard_cellchat=connected_full_CN[['Ligand', 'Receptor','significant']]
#print('Cellchat CN significant only')
#print(connected_full_CN_Jaccard_cellchat[connected_full_CN_Jaccard_cellchat['significant']==1].drop_duplicates())
connected_full_CN=connected_full_CN[['normal_receptor_list','cancer_ligand_list','significant']]
connected_full_CN['normal_receptor_list'] = connected_full_CN['normal_receptor_list'].apply(eval)
connected_full_CN['cancer_ligand_list'] = connected_full_CN['cancer_ligand_list'].apply(eval)
normal_receptor_expanded = connected_full_CN['normal_receptor_list'].apply(lambda x: pd.Series(x))
normal_receptor_expanded.columns = [f'normal_receptor_{i}' for i in normal_receptor_expanded.columns]
cancer_ligand_expanded = connected_full_CN['cancer_ligand_list'].apply(lambda x: pd.Series(x))
cancer_ligand_expanded.columns = [f'cancer_ligand_{i}' for i in cancer_ligand_expanded.columns]
df_expanded_CN = pd.concat([connected_full_CN.drop(['normal_receptor_list', 'cancer_ligand_list'], axis=1), normal_receptor_expanded, cancer_ligand_expanded], axis=1)

df_expanded_CN_cols = pd.DataFrame(columns=[f'column_{i}' for i in range(df_expanded_CN.shape[1])])
# Generate new column names
new_column_names = [str(i) for i in range(1, len(df_expanded_CN_cols.columns) + 1)]
# Assign new column names to the dataframe
df_expanded_CN.columns = new_column_names

#np.sum(df_expanded_CN['1']==df_expanded_NC['1'])

#print('Doing df_expanded_all=pd.concat((df_expanded_CN,df_expanded_NC)) for cellchat')

df_expanded_all=pd.concat((df_expanded_CN,df_expanded_NC))

#print('Dropping NAs from first column of df_expanded_all in cellchat')
df_expanded_all = df_expanded_all.dropna(subset=['1'])

features = df_expanded_all.drop('1', axis=1)
target = df_expanded_all['1']

# First, split the data into training (80%) and test (20%) sets
X_train, X_test, y_train, y_test = train_test_split(features, target, test_size=0.2, random_state=42, stratify=target)

# Print the shapes of the resulting DataFrames to verify the splits
#print(f'Training set: {X_train.shape}, {y_train.shape}')
#print(f'Test set: {X_test.shape}, {y_test.shape}')

# Initialize the Random Forest classifier
rf_model = RandomForestClassifier(n_estimators=100, random_state=42)

# Train the model on the training data
rf_model.fit(X_train, y_train)

# Make predictions on the test set
y_test_pred = rf_model.predict(X_test)

# Calculate accuracy
#test_accuracy = accuracy_score(y_test, y_test_pred)
#print(f'Cellchat test Accuracy: {test_accuracy:.4f}')

# Print classification report
print("Cellchat classification report:\n", classification_report(y_test, y_test_pred))

f1_cellchat = f1_score(y_test, y_test_pred, pos_label=1)

print('f1 score for CellChat is is '+ str(f1_cellchat))

# Print confusion matrix
#print("Cellchat confusion matrix:\n", confusion_matrix(y_test, y_test_pred))

print("Now doing random test")

####### Updated random test

num_iterations = 10
f1_scores = []
conf_matrices = []

print('Doing random iterations of '+str(num_iterations))
for i in range(num_iterations):
    print(f"\nIteration {i+1}:")
    # Generate a new random target vector
    random_target_vector = np.random.choice([0, 1], size=len(df_expanded_all['1']), p=[0.7, 0.3])
    # Split data into training (80%) and test (20%) sets
    X_train, X_test, y_train, y_test = train_test_split(features, random_target_vector, test_size=0.2, random_state=i, stratify=random_target_vector)
    # Initialize and train the model
    rf_model = RandomForestClassifier(n_estimators=100, random_state=i)
    rf_model.fit(X_train, y_train)
    # Make predictions
    y_test_pred = rf_model.predict(X_test)
    # Compute accuracy
    test_accuracy = accuracy_score(y_test, y_test_pred)
    #print(f"Random test Accuracy: {test_accuracy:.4f}")
    # Compute F1-score for class 1.0 and store it
    f1 = f1_score(y_test, y_test_pred, pos_label=1)
    f1_scores.append(f1)
    # Compute confusion matrix and store it
    conf_mat = confusion_matrix(y_test, y_test_pred)
    conf_matrices.append(conf_mat)
    # Print classification report & confusion matrix
    #print("Random classification report:\n", classification_report(y_test, y_test_pred))
    #print("Random confusion matrix:\n", conf_mat)

# Compute mean and standard deviation of F1 scores
f1_mean = np.mean(f1_scores)
f1_std = np.std(f1_scores)
# Compute average confusion matrix
avg_conf_matrix = np.mean(conf_matrices, axis=0).astype(int)
# Print final summary
print("\n====== FINAL SUMMARY ======")
print(f"Mean F1-score (class 1.0): {f1_mean:.4f}")
print(f"Standard Deviation of F1-score (class 1.0): {f1_std:.4f}")
print(f"Mean F1-score (class 1.0) + std: {f1_mean+f1_std:.4f}")
print(f"Mean F1-score (class 1.0) - std: {f1_mean-f1_std:.4f}")
#print("Average Confusion Matrix over 10 runs:\n", avg_conf_matrix)

###########################
######### LIANA results#########
###########################

print('Doing LIANA method - specificity')

######### Read in LIANA files (specificity and magnitude)

# Read the file as a single string and convert it into a list of tuples
if subset_type:
    with open(path + cell_subtype + '_Epithelial_'+ subset_type +'_specificity.txt', 'r') as f:
        data = ast.literal_eval(f.read().strip())  # Converts the string into a list of tuples
    print('Read in NC LIANA subset '+subset_type+' predictions')
else:
    with open(path + cell_subtype + '_Epithelial_specificity.txt', 'r') as f:
        data = ast.literal_eval(f.read().strip())  # Converts the string into a list of tuples
    print('Read in full NC LIANA specificity predictions')

# Convert to DataFrame with appropriate column names
NC_LIANA = pd.DataFrame(data, columns=['Ligand', 'Receptor'])

if subset_type:
# Do the same for the second file
    with open(path + 'Epithelial_' + cell_subtype + '_'+ subset_type +'_specificity.txt', 'r') as f:
        data = ast.literal_eval(f.read().strip())
    print('Read in CN LIANA specificity subset '+subset_type+' predictions')
else:
    with open(path + 'Epithelial_' + cell_subtype +'_specificity.txt', 'r') as f:
        data = ast.literal_eval(f.read().strip())
    print('Read in full CN LIANA specificity predictions')

CN_LIANA = pd.DataFrame(data, columns=['Ligand', 'Receptor'])

# Convert 'Ligand' and 'Receptor' columns to uppercase
NC_LIANA['Ligand'] = NC_LIANA['Ligand'].str.upper()
NC_LIANA['Receptor'] = NC_LIANA['Receptor'].str.upper()

CN_LIANA['Ligand'] = CN_LIANA['Ligand'].str.upper()
CN_LIANA['Receptor'] = CN_LIANA['Receptor'].str.upper()

connected_full_NC=pd.read_csv(path+cell_subtype+'_'+sample+'_concatenated_'+'NC_'+ct_level+'_'+version+'.csv')

print('Dropping NAs in LIANA results NC')
connected_full_NC = connected_full_NC.dropna(subset=['corr_connected'])

merged_LIANA_NC = connected_full_NC.merge(NC_LIANA, on=['Ligand', 'Receptor'], how='left', indicator=True)
connected_full_NC['significant'] = merged_LIANA_NC['_merge'].apply(lambda x: 1 if x == 'both' else 0)
#print('Merged LIANA_NC:')
#print(connected_full_NC[['Ligand', 'Receptor','significant']])
connected_full_NC_Jaccard_LIANA=connected_full_NC[['Ligand', 'Receptor','significant']]
#print('LIANA NC significant only')
#print(connected_full_NC_Jaccard_LIANA[connected_full_NC_Jaccard_LIANA['significant']==1].drop_duplicates())
connected_full_NC=connected_full_NC[['normal_ligand_list','cancer_receptor_list','significant']]
connected_full_NC['normal_ligand_list'] = connected_full_NC['normal_ligand_list'].apply(eval)
connected_full_NC['cancer_receptor_list'] = connected_full_NC['cancer_receptor_list'].apply(eval)
normal_ligand_expanded = connected_full_NC['normal_ligand_list'].apply(lambda x: pd.Series(x))
normal_ligand_expanded.columns = [f'normal_ligand_{i}' for i in normal_ligand_expanded.columns]
cancer_receptor_expanded = connected_full_NC['cancer_receptor_list'].apply(lambda x: pd.Series(x))
cancer_receptor_expanded.columns = [f'cancer_receptor_{i}' for i in cancer_receptor_expanded.columns]
df_expanded_NC = pd.concat([connected_full_NC.drop(['normal_ligand_list', 'cancer_receptor_list'], axis=1), normal_ligand_expanded, cancer_receptor_expanded], axis=1)

df_expanded_NC_cols = pd.DataFrame(columns=[f'column_{i}' for i in range(df_expanded_NC.shape[1])])
# Generate new column names
new_column_names = [str(i) for i in range(1, len(df_expanded_NC_cols.columns) + 1)]
# Assign new column names to the dataframe
df_expanded_NC.columns = new_column_names

#connected_full_CN=pd.read_csv(path+'Lymphocyte_ct_level0_raw_counts_v2_CN__main_known_ligrec_genes_part_1_complete_connected.csv')
#connected_full_CN=pd.read_csv(path+cell_subtype+'_'+ct_level+'_'+version+'_CN__main_known_ligrec_genes_part_1_complete_connected.csv')
connected_full_CN=pd.read_csv(path+cell_subtype+'_'+sample+'_concatenated_'+'CN_'+ct_level+'_'+version+'.csv')
#print('Dropping NAs in LIANA results')
connected_full_CN = connected_full_CN.dropna(subset=['corr_connected'])
merged_LIANA_CN = connected_full_CN.merge(CN_LIANA, on=['Ligand', 'Receptor'], how='left', indicator=True)
connected_full_CN['significant'] = merged_LIANA_CN['_merge'].apply(lambda x: 1 if x == 'both' else 0)
#print('Merged LIANA_CN:')
#print(connected_full_CN[['Ligand', 'Receptor','significant']])
connected_full_CN_Jaccard_LIANA_specificity=connected_full_CN[['Ligand', 'Receptor','significant']]
#print('LIANA CN significant only')
#print(connected_full_CN_Jaccard_LIANA[connected_full_CN_Jaccard_LIANA['significant']==1].drop_duplicates())
connected_full_CN=connected_full_CN[['normal_receptor_list','cancer_ligand_list','significant']]
connected_full_CN['normal_receptor_list'] = connected_full_CN['normal_receptor_list'].apply(eval)
connected_full_CN['cancer_ligand_list'] = connected_full_CN['cancer_ligand_list'].apply(eval)
normal_receptor_expanded = connected_full_CN['normal_receptor_list'].apply(lambda x: pd.Series(x))
normal_receptor_expanded.columns = [f'normal_receptor_{i}' for i in normal_receptor_expanded.columns]
cancer_ligand_expanded = connected_full_CN['cancer_ligand_list'].apply(lambda x: pd.Series(x))
cancer_ligand_expanded.columns = [f'cancer_ligand_{i}' for i in cancer_ligand_expanded.columns]
df_expanded_CN = pd.concat([connected_full_CN.drop(['normal_receptor_list', 'cancer_ligand_list'], axis=1), normal_receptor_expanded, cancer_ligand_expanded], axis=1)

df_expanded_CN_cols = pd.DataFrame(columns=[f'column_{i}' for i in range(df_expanded_CN.shape[1])])
# Generate new column names
new_column_names = [str(i) for i in range(1, len(df_expanded_CN_cols.columns) + 1)]
# Assign new column names to the dataframe
df_expanded_CN.columns = new_column_names

#np.sum(df_expanded_CN['1']==df_expanded_NC['1'])

print('Doing df_expanded_all=pd.concat((df_expanded_CN,df_expanded_NC)) for LIANA')

df_expanded_all=pd.concat((df_expanded_CN,df_expanded_NC))

#print('Dropping NAs from first column of df_expanded_all in LIANA')
df_expanded_all = df_expanded_all.dropna(subset=['1'])

features = df_expanded_all.drop('1', axis=1)
target = df_expanded_all['1']

# First, split the data into training (80%) and test (20%) sets
X_train, X_test, y_train, y_test = train_test_split(features, target, test_size=0.2, random_state=42, stratify=target)

# Initialize the Random Forest classifier
rf_model = RandomForestClassifier(n_estimators=100, random_state=42)

# Train the model on the training data
rf_model.fit(X_train, y_train)

# Make predictions on the test set
y_test_pred = rf_model.predict(X_test)

# Calculate accuracy
test_accuracy = accuracy_score(y_test, y_test_pred)
print(f'LIANA test Accuracy: {test_accuracy:.4f}')

# Print classification report
print("LIANA classification report (specificity):\n", classification_report(y_test, y_test_pred))

f1_liana_spec=f1_score(y_test, y_test_pred, pos_label=1)
print('f1 score for liana spec is '+ str(f1_liana_spec))

# Print confusion matrix
#print("LIANA confusion matrix (specificity):\n", confusion_matrix(y_test, y_test_pred))

print('Doing LIANA method - magnitude')

# Read the file as a single string and convert it into a list of tuples
if subset_type:
    with open(path + cell_subtype + '_Epithelial_'+ subset_type+'_magnitude.txt', 'r') as f:
        data = ast.literal_eval(f.read().strip())  # Converts the string into a list of tuples
    print('Read in NC LIANA magnitude subset '+subset_type+' predictions')
else:
    with open(path + cell_subtype + '_Epithelial_magnitude.txt', 'r') as f:
        data = ast.literal_eval(f.read().strip())  # Converts the string into a list of tuples
    print('Read in full LIANA magnitude NC LIANA predictions')


# Convert to DataFrame with appropriate column names
NC_LIANA = pd.DataFrame(data, columns=['Ligand', 'Receptor'])

# Do the same for the second file
if subset_type:
    with open(path + 'Epithelial_' + cell_subtype+ '_'+subset_type + '_magnitude.txt', 'r') as f:
        data = ast.literal_eval(f.read().strip())
    print('Read in CN LIANA magnitude subset '+subset_type+' predictions')
else:
    with open(path + 'Epithelial_' + cell_subtype+'_magnitude.txt', 'r') as f:
        data = ast.literal_eval(f.read().strip())  # Converts the string into a list of tuples
    print('Read in full LIANA magnitude CN LIANA predictions')

CN_LIANA = pd.DataFrame(data, columns=['Ligand', 'Receptor'])

# Convert 'Ligand' and 'Receptor' columns to uppercase
NC_LIANA['Ligand'] = NC_LIANA['Ligand'].str.upper()
NC_LIANA['Receptor'] = NC_LIANA['Receptor'].str.upper()

CN_LIANA['Ligand'] = CN_LIANA['Ligand'].str.upper()
CN_LIANA['Receptor'] = CN_LIANA['Receptor'].str.upper()

connected_full_NC=pd.read_csv(path+cell_subtype+'_'+sample+'_concatenated_'+'NC_'+ct_level+'_'+version+'.csv')

#print('Dropping NAs in LIANA results NC')
connected_full_NC = connected_full_NC.dropna(subset=['corr_connected'])

merged_LIANA_NC = connected_full_NC.merge(NC_LIANA, on=['Ligand', 'Receptor'], how='left', indicator=True)
connected_full_NC['significant'] = merged_LIANA_NC['_merge'].apply(lambda x: 1 if x == 'both' else 0)
#print('Merged LIANA_NC:')
#print(connected_full_NC[['Ligand', 'Receptor','significant']])
connected_full_NC_Jaccard_LIANA=connected_full_NC[['Ligand', 'Receptor','significant']]
#print('LIANA NC significant only')
#print(connected_full_NC_Jaccard_LIANA[connected_full_NC_Jaccard_LIANA['significant']==1].drop_duplicates())
connected_full_NC=connected_full_NC[['normal_ligand_list','cancer_receptor_list','significant']]
connected_full_NC['normal_ligand_list'] = connected_full_NC['normal_ligand_list'].apply(eval)
connected_full_NC['cancer_receptor_list'] = connected_full_NC['cancer_receptor_list'].apply(eval)
normal_ligand_expanded = connected_full_NC['normal_ligand_list'].apply(lambda x: pd.Series(x))
normal_ligand_expanded.columns = [f'normal_ligand_{i}' for i in normal_ligand_expanded.columns]
cancer_receptor_expanded = connected_full_NC['cancer_receptor_list'].apply(lambda x: pd.Series(x))
cancer_receptor_expanded.columns = [f'cancer_receptor_{i}' for i in cancer_receptor_expanded.columns]
df_expanded_NC = pd.concat([connected_full_NC.drop(['normal_ligand_list', 'cancer_receptor_list'], axis=1), normal_ligand_expanded, cancer_receptor_expanded], axis=1)

df_expanded_NC_cols = pd.DataFrame(columns=[f'column_{i}' for i in range(df_expanded_NC.shape[1])])
# Generate new column names
new_column_names = [str(i) for i in range(1, len(df_expanded_NC_cols.columns) + 1)]
# Assign new column names to the dataframe
df_expanded_NC.columns = new_column_names

#connected_full_CN=pd.read_csv(path+'Lymphocyte_ct_level0_raw_counts_v2_CN__main_known_ligrec_genes_part_1_complete_connected.csv')
#connected_full_CN=pd.read_csv(path+cell_subtype+'_'+ct_level+'_'+version+'_CN__main_known_ligrec_genes_part_1_complete_connected.csv')
connected_full_CN=pd.read_csv(path+cell_subtype+'_'+sample+'_concatenated_'+'CN_'+ct_level+'_'+version+'.csv')
#print('Dropping NAs in LIANA results')
connected_full_CN = connected_full_CN.dropna(subset=['corr_connected'])
merged_LIANA_CN = connected_full_CN.merge(CN_LIANA, on=['Ligand', 'Receptor'], how='left', indicator=True)
connected_full_CN['significant'] = merged_LIANA_CN['_merge'].apply(lambda x: 1 if x == 'both' else 0)
#print('Merged LIANA_CN:')
#print(connected_full_CN[['Ligand', 'Receptor','significant']])
connected_full_CN_Jaccard_LIANA_magnitude=connected_full_CN[['Ligand', 'Receptor','significant']]
#print('LIANA CN significant only')
#print(connected_full_CN_Jaccard_LIANA[connected_full_CN_Jaccard_LIANA['significant']==1].drop_duplicates())
connected_full_CN=connected_full_CN[['normal_receptor_list','cancer_ligand_list','significant']]
connected_full_CN['normal_receptor_list'] = connected_full_CN['normal_receptor_list'].apply(eval)
connected_full_CN['cancer_ligand_list'] = connected_full_CN['cancer_ligand_list'].apply(eval)
normal_receptor_expanded = connected_full_CN['normal_receptor_list'].apply(lambda x: pd.Series(x))
normal_receptor_expanded.columns = [f'normal_receptor_{i}' for i in normal_receptor_expanded.columns]
cancer_ligand_expanded = connected_full_CN['cancer_ligand_list'].apply(lambda x: pd.Series(x))
cancer_ligand_expanded.columns = [f'cancer_ligand_{i}' for i in cancer_ligand_expanded.columns]
df_expanded_CN = pd.concat([connected_full_CN.drop(['normal_receptor_list', 'cancer_ligand_list'], axis=1), normal_receptor_expanded, cancer_ligand_expanded], axis=1)

df_expanded_CN_cols = pd.DataFrame(columns=[f'column_{i}' for i in range(df_expanded_CN.shape[1])])
# Generate new column names
new_column_names = [str(i) for i in range(1, len(df_expanded_CN_cols.columns) + 1)]
# Assign new column names to the dataframe
df_expanded_CN.columns = new_column_names

#np.sum(df_expanded_CN['1']==df_expanded_NC['1'])

#print('Doing df_expanded_all=pd.concat((df_expanded_CN,df_expanded_NC)) for LIANA')

df_expanded_all=pd.concat((df_expanded_CN,df_expanded_NC))

#print('Dropping NAs from first column of df_expanded_all in LIANA')
df_expanded_all = df_expanded_all.dropna(subset=['1'])

# Filter out empty DataFrames before concatenation
#dfs_to_concat = [df for df in [df_expanded_CN, df_expanded_NC] if not df.empty]

# Concatenate only non-empty DataFrames
#df_expanded_all = pd.concat(dfs_to_concat, ignore_index=True) if dfs_to_concat else pd.DataFrame()

#print('df_expanded_all:\n')
#print(df_expanded_all)

#print('df_expanded_all')
#print(df_expanded_all)

features = df_expanded_all.drop('1', axis=1)
target = df_expanded_all['1']

# First, split the data into training (80%) and test (20%) sets
X_train, X_test, y_train, y_test = train_test_split(features, target, test_size=0.2, random_state=42, stratify=target)

# Initialize the Random Forest classifier
rf_model = RandomForestClassifier(n_estimators=100, random_state=42)

# Train the model on the training data
rf_model.fit(X_train, y_train)

# Make predictions on the test set
y_test_pred = rf_model.predict(X_test)

# Calculate accuracy
test_accuracy = accuracy_score(y_test, y_test_pred)
print(f'LIANA test Accuracy: {test_accuracy:.4f}')

# Print classification report
print("LIANA classification report (magnitude):\n", classification_report(y_test, y_test_pred))

# Print confusion matrix
#print("LIANA confusion matrix (magnitude):\n", confusion_matrix(y_test, y_test_pred))

f1_liana_mag=f1_score(y_test, y_test_pred, pos_label=1)
print('f1 score for liana mag is '+ str(f1_liana_mag))

# Extract patient identifier from the path
patient = os.path.basename(os.path.normpath(args.path))
    
# Example values from output (Replace these with actual computed values in your script)
#f1_home = 1.0  # Example from output
#f1_cellchat = 1.0  # Example from output
#f1_random = 0.3203  # Mean F1-score from random iterations
#f1_liana_spec = 1.0  # Example from specificity
#f1_liana_mag = 0.91  # Example from magnitude
    
#print(f"Mean F1-score (class 1.0) + std: {f1_mean+f1_std:.4f}")
#print(f"Mean F1-score (class 1.0) - std: {f1_mean-f1_std:.4f}")

f1_random_tuple = (round(f1_mean - f1_std, 4), round(f1_mean + f1_std, 4))  # Tuple format
    
if subset_type:
    print('Writing results for '+subset_type)
else:
    subset_type="Complete"
    
# Write results to file
write_results(patient, args.cell_subtype, subset_type, f1_home, f1_cellchat, f1_random_tuple, f1_liana_spec, f1_liana_mag)
print(f"Results appended for patient: {patient}, cell type: {args.cell_subtype}, subset type: {subset_type}")
    
# Extract Ligand-Receptor pairs for each method
liana_magnitude_pairs = set(
    zip(
        connected_full_CN_Jaccard_LIANA_magnitude.loc[connected_full_CN_Jaccard_LIANA_magnitude['significant'] == 1, 'Ligand'],
        connected_full_CN_Jaccard_LIANA_magnitude.loc[connected_full_CN_Jaccard_LIANA_magnitude['significant'] == 1, 'Receptor']
    )
)

cellchat_pairs = set(
    zip(
        connected_full_CN_Jaccard_cellchat.loc[connected_full_CN_Jaccard_cellchat['significant'] == 1, 'Ligand'],
        connected_full_CN_Jaccard_cellchat.loc[connected_full_CN_Jaccard_cellchat['significant'] == 1, 'Receptor']
    )
)

#home_pairs = set(
#    zip(
#        rows_of_interest['Ligand'],
#        rows_of_interest['Receptor']
#    )
#)

# Ensure rows_of_interest is defined before usage
if 'rows_of_interest' not in locals() and 'rows_of_interest' not in globals():
    rows_of_interest = pd.DataFrame()  # Initialize as an empty DataFrame

# Ensure home_pairs exists, even if rows_of_interest is empty or missing columns
if not rows_of_interest.empty and 'Ligand' in rows_of_interest and 'Receptor' in rows_of_interest:
    home_pairs = set(zip(rows_of_interest['Ligand'], rows_of_interest['Receptor']))
else:
    home_pairs = set()

# Ensure home_pairs exists, even if rows_of_interest is empty or missing columns
#if 'Ligand' in rows_of_interest and 'Receptor' in rows_of_interest and not rows_of_interest.empty:
#    home_pairs = set(zip(rows_of_interest['Ligand'], rows_of_interest['Receptor']))
#else:
#    home_pairs = set()

# Find common and unique pairs
common_pairs = liana_magnitude_pairs & cellchat_pairs & home_pairs
unique_to_liana_magnitude = liana_magnitude_pairs - (cellchat_pairs | home_pairs)
unique_to_cellchat = cellchat_pairs - (liana_magnitude_pairs | home_pairs)
unique_to_home = home_pairs - (liana_magnitude_pairs | cellchat_pairs)

# Compute Jaccard Index between methods
jaccard_liana_cellchat = jaccard_index(liana_magnitude_pairs, cellchat_pairs)
jaccard_liana_home = jaccard_index(liana_magnitude_pairs, home_pairs)
jaccard_cellchat_home = jaccard_index(cellchat_pairs, home_pairs)

# Create a DataFrame for output
comparison_df = pd.DataFrame({
    'Method': ['Liana', 'CellChat', 'Home'],
    'Total Pairs': [len(liana_magnitude_pairs), len(cellchat_pairs), len(home_pairs)],
    'Unique Pairs': [len(unique_to_liana_magnitude), len(unique_to_cellchat), len(unique_to_home)],
    'Common Pairs': [len(common_pairs), '', ''],
    'Jaccard Liana-CellChat': [jaccard_liana_cellchat, '', ''],
    'Jaccard Liana-Home': [jaccard_liana_home, '', ''],
    'Jaccard CellChat-Home': [jaccard_cellchat_home, '', '']
})

# Save all pairs in one document
all_pairs_df = pd.DataFrame(
    list(liana_magnitude_pairs | cellchat_pairs | home_pairs),
    columns=['Ligand', 'Receptor']
)
all_pairs_df['Method'] = all_pairs_df.apply(
    lambda row: ', '.join([
        'Liana_magnitude' if (row.Ligand, row.Receptor) in liana_magnitude_pairs else '',
        'CellChat' if (row.Ligand, row.Receptor) in cellchat_pairs else '',
        'Home' if (row.Ligand, row.Receptor) in home_pairs else ''
    ]), axis=1
)

# Save to a CSV file
if subset_type:
    comparison_df.to_csv(path+'method_comparison_'+subset_type+'.csv', index=False)
    all_pairs_df.to_csv(path+'all_ligand_receptor_pairs'+subset_type+'.csv', index=False)
else:
    comparison_df.to_csv(path+'method_comparison_.csv', index=False)
    all_pairs_df.to_csv(path+'all_ligand_receptor_pairs.csv', index=False)

# Print summary
print('Comparison between methods for '+cell_subtype+':\n')
print(comparison_df)

print('all pairs of each method for '+cell_subtype+':\n')
print(all_pairs_df)
print("\nSaved all ligand-receptor pairs in 'all_ligand_receptor_pairs.csv'.")

#########
####THE END####
#########




