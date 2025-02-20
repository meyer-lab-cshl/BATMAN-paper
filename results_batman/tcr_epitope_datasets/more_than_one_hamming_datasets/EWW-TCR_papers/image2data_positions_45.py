import numpy as np
import pandas as pd
from PIL import Image
import itertools
from matplotlib import pyplot as plt

'''Extract RGB values'''
# Load image
im = Image.open('all_heatmaps.jpg','r') 
pix = list(im.getdata())
width, height = im.size
pixel_values = np.array(pix).reshape((width, height, 3))

# Set up a grid to extract RGB values for each box in heatmap

# For mutations in positions 4 and 5
xmin = 203.5
xmax = 851.5
ymin = 162.5
ymax = 824.5
xpts = 11
ypts = 12

grid = list(itertools.product(np.linspace(xmin,xmax,xpts),
                              np.linspace(ymin,ymax,ypts)))

# Extract RGB values
rgb = np.zeros((xpts*ypts,3))

for pixel in np.arange(xpts*ypts):
    rgb[pixel,:] = im.getpixel(grid[pixel])

# Reshape to heatmap size y channel
rgb = np.reshape(rgb, (xpts,ypts,3))

'''Colorbar'''
# sample equally spaced 100 collinear points between upper (6) & lower (-6) ticks

xmin = 3502
xmax = 3503
ymin = 2675
ymax = 3162.5
xpts_cb= 1
ypts_cb = 100

grid = list(itertools.product(np.linspace(xmin,xmax,xpts_cb),
                              np.linspace(ymin,ymax,ypts_cb)))

# Extract RGB values
cb_rgb = np.zeros((ypts_cb,3))

for pixel in np.arange(ypts_cb):
    cb_rgb[pixel,:] = im.getpixel(grid[pixel])

# Interpolate logFC scale
rgb2log_FC = np.linspace(6,-6,ypts_cb)

# RGB value of peptides not screened
rgb_unscreened = im.getpixel((3807,3135))

'''Assign logFC values based on RGB'''
mutant_logFC = pd.DataFrame(columns=['peptide','peptide_activity'])

# logFC normalization from min and max of 1 Hamming scan of EWW peptide
logFC_min = 0
logFC_max = 5.2

#mutation positions
Px = 4
Py = 5

# mutated AAs
index_peptide = "EWWRSGGFSF"
aa_order_x = ['A','C','E','F','G','I','K','M','N','P','T']
aa_order_y = ['T','R','P','N','M','K','I','G','F','E','C','A']

for aa_x in np.arange(xpts):
    for aa_y in np.arange(ypts):         
        # Get RGB values
        seq_rgb = rgb[aa_x,aa_y,:] #Read recorded rgb value
        # Find closest RGB combination in the colorbar
        cb_loc=np.argmin(np.sum(abs(cb_rgb-np.tile(seq_rgb,(len(cb_rgb),1))),axis=1))
        seq_logFC = rgb2log_FC[cb_loc] #assign logFC
        # Normalize logFC
        seq_logFC = max(0.0,(seq_logFC-logFC_min)/(logFC_max-logFC_min))
        
        # Detect unscreened
        if np.sum(abs(seq_rgb-rgb_unscreened))<np.min(np.sum(abs(cb_rgb-np.tile(seq_rgb,(len(cb_rgb),1))),axis=1)):
           seq_logFC = np.nan 
           continue
        # Mutant seq
        peptide = list(index_peptide)
        peptide[Px-1] = aa_order_x[aa_x]
        peptide[Py-1] = aa_order_y[aa_y]
        peptide = ''.join(peptide)
        
        # Save data
        data_add = pd.DataFrame({'peptide':[peptide],
                                 'peptide_activity':[seq_logFC]})
        mutant_logFC = pd.concat([mutant_logFC,data_add])
        

# Export data
mutant_logFC.to_excel('data_from_heatmap_45.xlsx',index=None)
       
        





















