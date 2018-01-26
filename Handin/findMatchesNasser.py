#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 23:31:13 2017

@author: nannabarnkob
"""


import pandas
import os 
import numpy as np 
import sys
import pdb
import pylab as plt
from matplotlib_venn import venn3, venn3_circles
import re

def findMatchesNasser(hallmarkFile, variantFile):
    df_HM = pandas.read_excel(hallmarkFile)
    with open(variantFile) as f:
       content = f.readlines()
       # you may also want to remove whitespace characters like `\n` at the end of each line
       geneSymbols_variants = [x.strip() for x in content] 

    geneSymbols_variants=[x.strip('"') for x in geneSymbols_variants[0].split(',')]
#    print(geneSymbols_variants)
    
    #get the values for gene symbol columns 
    hallmarks = list(df_HM.columns)
    HM_dct = {}
    match_dct = {}
    HMsizes = []
    labels = []
    unreportedHMs = []
    geneMatches= []
    title = re.search('(.*)_geneList',str(variantFile))
    outfile = open(title.group(1)+'_HMmatches_Nasser.txt','w')
    outfile.writelines('Hallmark \t Count \t Gene symbols \n')

    for hallmark in hallmarks:
 #       print('working on', hallmark)
        HM_dct = [str(x) for x in df_HM[hallmark].values]
        match= list(set(HM_dct).intersection(geneSymbols_variants))
        #pdb.set_trace()
        #match = set(list(df_HM[hallmark])).intersection(geneSymbols_variants)
        outfile.writelines([hallmark, '\t',str(len(match)),'\t', ','.join(map(str, list(match))), '\n'])
        [geneMatches.append(x) for x in match]
        if len(match) != 0:
            HMsizes.append(len(match))
            labels.append(hallmark)
        else:
            unreportedHMs.append(hallmark)
    
    # count statistics
    outfile.writelines(['Total input genes \t', str(len(geneSymbols_variants)),'\n'])
    outfile.writelines(['Unique genes mapped to Hallmarks \t', str(len(set(geneMatches))),'\n'])
    outfile.writelines(['Number of gene observations \t', str(np.sum(HMsizes)),'\n'])

hallmarkFile='/home/projects/dp_00005/data/nanbar/AnnotationScripts/Hallmarks_of_Cancer_TGen.xlsx'

for fileName in sys.argv[1:]:
    variantFile=fileName
    findMatchesNasser(hallmarkFile, variantFile)
