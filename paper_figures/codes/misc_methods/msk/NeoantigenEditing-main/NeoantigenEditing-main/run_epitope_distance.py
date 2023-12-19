import pandas as pd
import numpy
from scipy.io import savemat
from EpitopeDistance import EpitopeDistance

# Load epitope list
epitope_seqs = pd.read_excel('TCR_epitope_database.xlsx')['peptide'];
epitope_seqs=epitope_seqs.to_numpy();
epitope_seqs=numpy.unique(epitope_seqs);

# Filter 9-mers
s=pd.Series(epitope_seqs);
seqs=s.loc[s.str.len() == 9].to_numpy()

# Find pairwise epitope-distances from Immunoediting paper

d=numpy.zeros((len(seqs), len(seqs)));

epidist = EpitopeDistance()

for seq1 in numpy.arange(0,len(seqs),1):
    print(seq1)
    for seq2 in numpy.arange(0,len(seqs),1):

        d[seq1,seq2]=epidist.epitope_dist(seqs[seq1], seqs[seq2])

seq_dist=pd.DataFrame(d)
seq_dist.index=seqs;
seq_dist.columns=seqs;
seq_dist.to_csv('distances.csv')
#Remember to add the column name "seqs" in the csv upper left corner
