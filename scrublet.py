#see https://github.com/AllonKleinLab/scrublet for complete documentation on this tool

#Install 'Scrublet'
pip install scrublet

#Start Python (3.8)
python

#############
# In Python #
#############

#Load required packages
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os

#Set file directory
input_dir = 'input_dir' #set the correct file directory

#Load sparse matrix from R
counts_matrix = scipy.io.mmread(input_dir + '/Scrublet_expression_sparse_t.mtx')

#Run scrublet with default parameters
scrub = scr.Scrublet(counts_matrix) #expected_doublet_rate=0.10

#Extract scores and predictions from Scrublet output with default parameters
doublet_scores, predicted_doublets = scrub.scrub_doublets() #min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30

#Write out both files as .txt
np.savetxt(input_dir + '/Doublet_scores.txt', doublet_scores, delimiter="\t", fmt="%0.4f")
np.savetxt(input_dir + '/Predicted_doublets.txt', predicted_doublets, delimiter="\t", fmt="%d")

#histogram visualization

#Examine histogram and save in project folder
scrub.plot_histogram()[0].savefig(input_dir + '/test.png')
