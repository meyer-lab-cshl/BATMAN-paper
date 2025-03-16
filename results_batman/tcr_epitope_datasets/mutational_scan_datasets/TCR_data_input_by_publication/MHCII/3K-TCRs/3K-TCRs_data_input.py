import numpy as np
import pandas as pd
from PIL import Image
import itertools
from matplotlib import pyplot as plt

# Fixed data
tcr_name_list = ['B3K2.8',
                 'B3K503','B3K506','B3K525','B3K508','B3K5034','B3K5038',
                 'TCR72-18', 'TCR72-27', 'TCR73-37', 'TCR73-130',
                 'TCR75-23', 
                 'TCR77-55','YAe5-100.5','YAe5-62.8',
                 'YAe5-15.8',
                 'TCR75-1','TCR74-1','YAe9-14.1', 'TCR75-49', 'TCR2W1S4-43',
                 'YAe9-99','YAe6-13.2','TCR74-34','TCR75-21','TCR75-54',
                 'YAe10-106.3','TCR75-10','TCR74-36','TCR74-18','TCR74-22','TCR2W1S12-20.4',
                 'TCR2W1S12-146.3', 'YAe5-97', 'TCR74-24','TCR74-5','TCR75-13','TCR75-7','TCR2W1S12-118.6',
                 'TCR74-27','YAe10-20.3',
                 'TCR3K-36']

cdr3a_list = ['',
              '','CALVISNTNKVVFG','','CAASKGNNRIFFG','','',
              '','','','',
              '',
              '','','CAANSGTYQRFG',
              '',
              'CAASKNAPRFG','','','','',
              '','','','','',
              '','','','','','CAASDNRIFFG',
              '','','','','','','',
              '','',
              '']

cdr3b_list = ['',
              '','CASIDSSGNTLYFG','','CAWSLGGVAETLYFG','','',
              '','','','',
              '',
              '','','CASGDFWGDTLYFG',
              '',
              'CASSDSSSYEQYFG','','','','',
              '','','','','',
              '','','','','','CASGDAWGYEQYFG',
              '','','','','','','',
              '','',
              '']

trav_list = ['',
              '','4.1','','2.3','','',
              '','','','',
              '',
              '','','4S6/12',
              '',
              '2S3','4.11','','','',
              '','','','','',
              '','','','','','2S3',
              '','','','','','','',
              '','',
              '3.4']

traj_list = ['',
              '','27 (34)','','24 (31)','','',
              '','','','',
              '',
              '','','11',
              '',
              '35 (43)','','','','',
              '','','','','',
              '','','','','','24',
              '','','','','','','',
              '','',
              '']

trbv_list = ['',
              '','8.1','','14','','',
              '','','','',
              '',
              '','','8.2',
              '',
              '8.1','14','','','',
              '','','','','',
              '','','','','','8.2',
              '','','','','','','',
              '','',
              '8.2']

trbj_list = ['',
              '','2.3','','2.3','','',
              '','','','',
              '',
              '','','2.4',
              '',
              '2.6','','','','',
              '','','','','',
              '','','','','','2.6',
              '','','','','','','',
              '','',
              '']

'''Extract RGB values'''
# Load image
im = Image.open('heatmaps_all_tcr.jpg','r') 
pix = list(im.getdata())
width, height = im.size
pixel_values = np.array(pix).reshape((width, height, 3))

# Set up a grid to extract RGB values for each box in heatmap
xpts = np.concatenate((np.array([392]),
                       np.linspace(454,656,6),
                       np.linspace(798,919,4),
                       np.array([1192]),
                       np.linspace(1253,1334,3),
                       np.array([1397]),
                       np.linspace(1458.5,1618.5,5),
                       np.linspace(1680,2082,11),
                       np.linspace(2144,2385,7),
                       np.array([2447,2487,2713])
                       ))

ypts = np.concatenate((np.linspace(610.5,1245,19),
                       np.linspace(610.5,1245,19)+706,
                       np.linspace(610.5,1245,19)+2*706,
                       np.linspace(610.5,1245,19)+3*706,
                       np.linspace(610.5,1245,19)+4*706
                       ))

grid = list(itertools.product(xpts,
                              ypts))

# Extract RGB values
rgb = np.zeros((len(xpts)*len(ypts),3))

for pixel in np.arange(len(xpts)*len(ypts)):
    rgb[pixel,:] = im.getpixel(grid[pixel])

# Reshape to heatmap size y channel
rgb = np.reshape(rgb, (len(xpts),len(ypts),3))

rgb_classes = np.stack((np.array([135,45,18]),# red
                              np.array([238,228,42]),# yellow
                              np.array([255,255,255]),# white
                              ))
activity_values = [0.02,0.25,0.75]

aa_list = ['G','P','S','T','A','V','I','L','M','C','F','Y','W','N','Q','D','K','R','H',
           'G','P','S','T','A','V','I','L','M','C','F','Y','W','N','D','E','K','R','H',
           'G','P','S','T','A','V','I','L','M','C','F','Y','W','N','Q','D','E','R','H',
           'G','P','S','T','A','V','I','L','M','C','F','Y','W','N','Q','D','E','R','H',
           'G','P','S','T','A','V','I','L','M','C','F','Y','W','N','Q','D','E','R','H']

mutation_position_list = np.concatenate((1*np.ones(19),
                                         3*np.ones(19),
                                         4*np.ones(19),
                                         6*np.ones(19),
                                         9*np.ones(19)
                                         ))

index_peptide = "FEAQKAKANKAVD"
# Initialize df
all_data=pd.DataFrame()

for tcr in np.arange(len(xpts)):
    for mutant_aa in np.arange(len(ypts)):
        
        #Skip "ND" values
        if tcr_name_list[tcr] in ['B3K5034','TCR73-130'] and int(mutation_position_list[mutant_aa])==1 and aa_list[mutant_aa]!='A':
            continue
            
        
        # Record TCR, MHC, organism etc info
        tcr_data = pd.DataFrame({
        'tcr_name':[tcr_name_list[tcr]],
        'va':["NA"],
        'vb':["NA"],
        'cdr3a':[cdr3a_list[tcr]],
        'cdr3b':[cdr3b_list[tcr]],
        'trav':[trav_list[tcr]],
        'traj':[traj_list[tcr]],
        'trbv':[trbv_list[tcr]],
        'trbd':["NA"],
        'trbj':[trbj_list[tcr]],
        'assay':["T cell proliferation"],
        'tcr_source_organism':["mouse"],
        'index_peptide':["FEAQKAKANKAVD"],
        'mhc':["I-Ab"],
        'pmid':["16051149"],
        'peptide_type':["unspecified"]
        }) 
        
        # Get RGB values
        seq_rgb = rgb[tcr,mutant_aa,:] #Read recorded rgb value
        # Find closest RGB combination in the colorbar and assign value
        peptide_activity=activity_values[
        np.argmin(np.sum(abs(rgb_classes-np.tile(seq_rgb,(len(rgb_classes),1))),axis=1))]
        
        # Construct mutant sequence
        peptide = list(index_peptide)
        peptide[int(mutation_position_list[mutant_aa])] = aa_list[mutant_aa]
        peptide = ''.join(peptide)        
            
        # Add data
        data_add = pd.DataFrame({'peptide':[peptide],
                                 'peptide_activity':[peptide_activity]})
        data_add = pd.concat([tcr_data,data_add],axis=1)        
        all_data = pd.concat([all_data,data_add])
        
    # Add unmutated peptide
    data_add = pd.DataFrame({'peptide':[index_peptide],
                             'peptide_activity':[1]})
    data_add = pd.concat([tcr_data,data_add],axis=1)        
    all_data = pd.concat([all_data,data_add])
    
        
all_data=all_data.drop_duplicates().reset_index(drop=True)
# Export data
all_data.to_excel('epitope_data_3K-TCR.xlsx', index=False)

        
        



















