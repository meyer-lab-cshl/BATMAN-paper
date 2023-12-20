# BATMAN: Bayesian inference of Activation of TCR by Mutant Antigens
A single T Cell Receptor (TCR) can recognize a diverse variety of peptides, an essential property known as TCR cross-reactivity. Predicting which peptides a TCR cross-reacts to is critical for numerous applications, including predicting viral escape, cancer neoantigen immunogenicity, autoimmunity, and off-target toxicity of T-cell-based therapies. But predicting TCR activation is challenging due to the lack of both unbiased benchmarking datasets and computational methods that are sensitive to small mutations to an epitope. To address these challenges, we curated a comprehensive database encompassing complete single-amino-acid mutational assays of 10,750 TCR-peptide pairs, centered around 14 immunogenic epitopes against 66 TCRs. We then developed an interpretable Bayesian model, called BATMAN, that can predict the set of epitopes that activate a TCR. When validated on our database, BATMAN outperforms existing methods by 20% and reveals important biochemical predictors of TCR-peptide interactions.

BATMAN predicts TCR activation by mutant peptides based on their distances to the TCR's index peptide. The peptide-to-index distance is a product of a learned positional weight profile vector, corresponding to effects of mutated residues at different positions in the sequence, and a learned AA substitution distance from the index peptide amino acid to the mutant amino acid.

![BATMAN schematic diagram showing that it integrates mutational scan datasets across many TCRs to build a hierarchical Bayesian inference model. BATMAN infers
hyperparameters from the training database and uses them to generate prior distributions for cross-TCR AA distance and TCR-specific
positional weights, which are multiplied and used as a predictor of TCR activation by a given mutant](BATMAN_schematic.jpg)

BATMAN can be trained in two modes: (1) within-TCR, where the train and test peptides are associated with the same TCR, and BATMAN-inferred positional weight profiles are TCR-specific, and (2) leave-one-TCR-out, where peptides are tested for activation of a TCR left out of the training data, and BATMAN-inferred positional weight profile is common across all TCRs.

For more information, refer to our preprint! For an interactive tutorial and test input, refer to our [jupyter notebook](https://github.com/meyer-lab-cshl/BATMAN-paper/blob/main/run_batman/pyBATMAN_Tutorial.ipynb).

# paper_figures
Folder contains all codes and raw data to reproduce figures in the BATMAN paper. The excel file paper_figures/data/TCR_epitope_database.xlsx is the main file containing all the TCR-pMHC mutational scan data we curated.

# pybatman
Folder contains codes for the Python package pyBATMAN, hosted at https://pypi.org/project/pybatman/.

# run_batman
Folder contains a test input data file and a python test script that trains BATMAN on the input data file and plots results. Please pip install pyBATMAN, preferably in a new conda environment with python=3.11, before running the script. The folder also has an interactive Jupyter notebook that trains and validates pyBATMAN on the test data and visualizes the results.
