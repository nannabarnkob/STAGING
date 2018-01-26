#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 12:56:35 2018

@author: nannabarnkob
This script maps gene variants to iallmark capabilities based on the dataset by Knijnenburg et al. 
available through Synapse.org with the identifier syn4216888 (http://dx.doi.org/10.7303/syn4216888)  
A pie chart is also produced. 
"""
# for everything 
import pandas as pd
import os
import numpy as np
import re
# for handling command-line input 
import sys
# for plotting 
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

def findMatchesKnijnenburg(hallmarkFile,inputVariants):
    df = pd.read_table(hallmarkFile, index_col=0)
    cols = df.columns
    print(inputVariants)
    sampleName = re.search('(.*)_geneList',str(inputVariants))
    # reading input file to python list 
    with open(inputVariants) as f:
       content = f.readlines()
       content = content[0].split(',')
       content = [x.replace('\"','') for x in content]
       # making sure to remove whitespace characters like `\n` at the end of each line
       variantGenes = [x.strip() for x in content]
    
    # create set of input variants for comparison   
    variantGenesSet = set(variantGenes)

    # writes to file 
    outfile = open(sampleName.group(1)+'_'+'HMmatches_Knijnenburg.txt','w')

    # specify header 
    outfile.writelines('Hallmark \t Count \t Gene symbols \n')
    HMsizes = []
    labels = []
    unreportedHMs = [] 
    geneMatches= []
    
    for hallmark in cols:
        ts = pd.Series(df[hallmark])
        HMgenes = list(ts.index[ts.nonzero()])
        match = list(set(HMgenes).intersection(variantGenesSet))
        outfile.writelines([hallmark, '\t', str(len(match)),'\t', ','.join(map(str, list(match))),'\n'])
        [geneMatches.append(x) for x in list(set(HMgenes).intersection(variantGenes))]
        # specify values to plot 
        if len(match) != 0:
            HMsizes.append(len(match))
            labels.append(hallmark)
        else: 
            unreportedHMs.append(hallmark)
    # plot settings 
    sizes = HMsizes

    # count statistics
    outfile.writelines(['Total input genes \t', str(len(variantGenes)),'\n'])
    outfile.writelines(['Unique genes mapped to Hallmarks \t', str(len(set(geneMatches))),'\n'])
    outfile.writelines(['Number of gene observations \t', str(np.sum(sizes)),'\n'])    

#    # PLOTTING BAR CHART 
#    index = np.arange(len(labels))
#    fig = plt.figure()
#    
#    # define position of labels 
#    ha = ['right', 'center', 'left']
#    # create subplot for axis 
#    ax = fig.add_subplot(111)
#    ax.set_xticks(index)
#    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
#    for n in range(3):
#        ax.set_xticklabels(labels, ha=ha[n])
#    for tick in ax.get_xticklabels():
#        tick.set_rotation(310)
#    # create body of plot 
#    plt.bar(index, sizes, width=0.7, align="center")
#    plt.title(sampleName.group(1))
#    # save barchart figure to file
#    plt.savefig(inputVariants+'_BarPlot.png')
#
#    # PLOTTING OF PIE CHART 
#    # make one part stand out (this will currently be the first HM per default)
#    explode = (0.1,) + (0,)*(len(HMsizes)-1)
#    
#    # specify colormap and calculate colors 
#    number = len(HMsizes)
#    cmap = plt.get_cmap('RdBu')
#    colors = [cmap(i) for i in np.linspace(0, 1, number)]
#
#    # set size 
#    plt.figure(figsize=(10,10))
#    plt.axis('equal')
#    
#    patches, texts, autotexts = plt.pie(sizes, explode=explode, labels=labels, colors=colors,
#            autopct='%1.1f%%', shadow=True, startangle=140)
#    # add legend 
#    plt.legend(patches, labels, loc="best")
#    
#    [i.set_fontsize(15) for i in texts]
#    [i.set_fontsize(20) for i in autotexts]
#    autotexts[-1].set_color('white')    
#    
#    plt.title(sampleName.group(1)) 
#    # save figure to file 
#    plt.savefig(inputVariants+'_PieChart.png')

hallmarkFile='/home/projects/dp_00005/data/nanbar/AnnotationScripts/Genes2Hallmarks.tsv'

for fileName in sys.argv[1:]:
      inputVariants=fileName
      findMatchesKnijnenburg(hallmarkFile, inputVariants)
