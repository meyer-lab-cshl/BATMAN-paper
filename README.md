# BATMAN: Bayesian inference of Activation of TCR by Mutant Antigens
A single T Cell Receptor (TCR) can recognize a diverse variety of peptides, an essential property known as TCR cross-reactivity. Predicting which peptides a TCR cross-reacts to is critical for numerous applications, including predicting viral escape, cancer neoantigen immunogenicity, autoimmunity, and off-target toxicity of T-cell-based therapies. But predicting TCR activation is challenging due to the lack of both unbiased benchmarking datasets and computational methods that are sensitive to small mutations to an epitope. To address these challenges, we curated a comprehensive database encompassing complete single-amino-acid mutational assays of 10,750 TCR-peptide pairs, centered around 14 immunogenic epitopes against 66 TCRs. We then developed an interpretable Bayesian model, called BATMAN, that can predict the set of epitopes that activate a TCR. When validated on our database, BATMAN outperforms existing methods by 20% and reveals important biochemical predictors of TCR-peptide interactions.

# paper_figures
Folder contains all codes and raw data to reproduce figures in the BATMAN paper. The excel file paper_figures/data/TCR_epitope_database.xlsx is the main file containing all the TCR-pMHC mutational scan data we curated.

# pybatman
Folder contains codes for the Python package pyBATMAN, hosted at https://pypi.org/project/pybatman/.

# run_batman
Folder contains a test input data file and a python test script that trains BATMAN on the input data file and plots results. Please pip install pyBATMAN, preferably in a new conda environment with python=3.11, before running the script.
