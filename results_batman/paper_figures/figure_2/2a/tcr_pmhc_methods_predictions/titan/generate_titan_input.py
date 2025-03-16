import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# Import BATMAN database
tcr_pmhc_data = pd.read_excel(''.join(['../../../../../tcr_epitope_datasets/',
          'mutational_scan_datasets/train_test_data_folds.xlsx']),index_col=0)

# Subset selected TCRs
TCR_list_9mer = ['a3a','TIL1383I','TCR81-14','18A2','NYE-S1','TCR6','TCR7',
                 'A6','T1','FLT3DY','A23','TCR2-T','R24','868Z11'] 

TCR_list_10mer = ['TCR-1E6', 'APN-TCR', 'EWW-TCR', 'A11Va']

TCR_list_mhcii = ['TCR-F5','TCR-3598-2','MBP-TCR','B3K508']

TCR_list = TCR_list_9mer + TCR_list_10mer + TCR_list_mhcii

tcr_pmhc_data = tcr_pmhc_data[tcr_pmhc_data.tcr.isin(TCR_list)]

# Subset to TCRs with full TCRb seq
tcr_pmhc_data=tcr_pmhc_data[pd.isna(tcr_pmhc_data.trbv)==False]
tcr_pmhc_data=tcr_pmhc_data[pd.isna(tcr_pmhc_data.trbj)==False]
tcr_pmhc_data=tcr_pmhc_data[pd.isna(tcr_pmhc_data.cdr3b)==False]

# Export to generate TCR sequence list (used: https://tcrmodel.ibbr.umd.edu/)
# and IMGT https://www.imgt.org/IMGTrepertoire/Proteins/alleles/index.php?species=Homo%20sapiens&group=TRBV&gene=TRBV10-3
# To have the first 3 TRBV AAs to conform to the convention of https://github.com/PaccMann/TITAN/blob/main/tutorial/data/tcr_full.csv
#tcr_pmhc_data[['tcr','trbv','trbj','cdr3b','tcr_source_organism'
#               ]].drop_duplicates().to_excel('TCR_data.xlsx',index=None)

# Assign TCR and epitope IDs
tcr_pmhc_data['tcr_id'] = tcr_pmhc_data.groupby('tcr').ngroup() +3000 
#3000 added so that TCR and epitope IDs do not overlap

tcr_pmhc_data['epitope_id'] = tcr_pmhc_data.groupby('peptide').ngroup() + 3

# Load and export TCR ID file
tcr_ids = tcr_pmhc_data[['tcr','tcr_id']].drop_duplicates()
tcr_seqs = pd.read_excel('TCR_data.xlsx',index_col=0)
tcr_ids['tcr_seq'] = tcr_seqs.TCRb_full[tcr_ids['tcr']].to_numpy()
tcr_ids[['tcr_seq','tcr_id']].to_csv('tcr_seq_id.csv',sep='\t',index=None,header=None)

# Export epitope ID file to smi
epitope_ids = tcr_pmhc_data[['peptide','epitope_id']].drop_duplicates()

with open('epitope_seq_id.smi', 'w') as f:
    for i in np.arange(len(epitope_ids)):
        ID = epitope_ids.iloc[i,1]
        
        # Convert AA seq to smiles
        seq = epitope_ids.iloc[i,0]
        mol = Chem.MolFromSequence(seq)
        smi = Chem.MolToSmiles(mol)
        
        # Add line to file
        f.write(f"{smi}\t{ID}\n")  # Replace 'Name' with the appropriate column name

epitope_ids.to_csv('epitope_seq_id.csv',index=None,sep='\t',header=None)

# Export test file
# Format data according to TITAN format
#Example: https://github.com/PaccMann/TITAN/tree/main/tutorial/data/test_small.csv
tcr_pmhc_data = tcr_pmhc_data[['epitope_id','tcr_id']]

tcr_pmhc_data.columns = ['ligand_name','sequence_id']
tcr_pmhc_data['label'] = 0 #dummy
tcr_pmhc_data['dummy_id']=tcr_pmhc_data['ligand_name'].to_numpy()
tcr_pmhc_data = tcr_pmhc_data[['dummy_id','ligand_name','sequence_id','label']]
tcr_pmhc_data.columns=['','ligand_name','sequence_id','label']


# Export data
tcr_pmhc_data.to_csv('titan_input_peptides.csv',index=None)
