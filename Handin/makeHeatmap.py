#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 15:26:00 2018

@author: nannabarnkob
"""
# resource for dataframe handling 
import pandas as pd 
# for regular expressions 
import re 
# data manipulation / normalizatoin 
from sklearn.preprocessing import MinMaxScaler
# for handling command-line input
import sys 
# resources for plotting 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import types
import pdb

def makeHeatmap(germlineVariants, somaticVariants, subtractedVariants, mutectVariants, combinedCounts):
    # find sample we're working with
    title = re.search('(.*)_normal_HMmatches_combined.txt',str(germlineVariants))
    pretty_title = title.group(1).replace('_',' all ')
    pretty_title = pretty_title.replace('DelDam', 'intersect')
    print('Statistics for '+pretty_title)
    # reading and preparing data
    df_germline = pd.read_table(germlineVariants, index_col=0)
    df_somatic = pd.read_table(somaticVariants, index_col=0)
    df_subtracted =  pd.read_table(subtractedVariants, index_col=0)
    df_mutect = pd.read_table(mutectVariants, index_col=0)
    df_totalCounts = pd.read_table(combinedCounts, index_col=0)
    # All counts are 'extracted' 
    germline_count = df_germline.loc[:,' Count ']
    somatic_count = df_somatic.loc[:,' Count ']
    subtracted_count = df_subtracted.loc[:,' Count ']
    mutect_count = df_mutect.loc[:,' Count ']
    #print(df_totalCounts)
    count_series = [germline_count, somatic_count, subtracted_count, mutect_count]
    for count in count_series:
        if sum(count) == 0:
            print("No somatic variants found, cannot construct heatmap")
    #series = [germline_count, somatic_count]
    data = pd.concat(count_series, axis=1)
    data.columns=['Germline \n ['+str(germline_count.iloc[0:10].sum(axis=0))+' / '+str(germline_count.loc['Total unique genes mapped to hallmarks '])+']', 'Tumor \n ['+str(somatic_count.iloc[0:10].sum(axis=0))+' / '+str(somatic_count.loc['Total unique genes mapped to hallmarks '])+']',
        'Subtracted \n ['+str(subtracted_count.iloc[0:10].sum(axis=0))+' / '+str(subtracted_count.loc['Total unique genes mapped to hallmarks '])+']', 'Mutect \n ['+str(mutect_count.iloc[0:10].sum(axis=0))+' / '+str(mutect_count.loc['Total unique genes mapped to hallmarks '])+']']
    #print(data.index)
    # when using the combined dataset, it is not neccessary to drop 'Other' (only for Hallmark tables made using K-mapping data) 
    data = data.drop(['Total unique genes mapped to hallmarks ', 'Number of gene mappings '])    
    print("Genes pr. hallmark before any normalization:")
    print(data)
    # normalization with respect ot number of genes in each category 
    data = data.divide(df_totalCounts.loc[:,' Count '], axis='index')
    print("Genes pr. hallmark  after total number of genes pr. hallmark normalization:")
    print(data)
    # normalization with respect to number of genes in each analysis 
    # old normalization based on rescaling ( x' = x - x[min] / x[max] - x[min] ) 
    # scaler = MinMaxScaler()
    # scaled_values = scaler.fit_transform(data) 
    # data.loc[:,:] = scaled_values  
    # Now substituted for normalization with respect to sum of each columns (analysis) 
    data = data.divide(data.sum(axis=0),axis='columns')
    print("Genes pr. hallmark after normalization espect to number of genes in each analysis:")
    print(data)
    # Finally plot in heatmap  
    plt.figure(figsize=(10,10))
    sns.set(font_scale=1.8, font='serif') 
    r = sns.heatmap(data, cmap='BuPu', square=True)
    ha = ['left', 'center', 'right']
    #for n in range(3):
    #    r.set_xticklabels(data.columns,ha=ha[n])
    for tick in r.get_xticklabels():
        tick.set_rotation(90)
    r.set_title(pretty_title+' variants')
    plt.gcf().subplots_adjust(left=0.45)
    plt.xticks(rotation=90)
    plt.savefig(title.group(1)+'_heatMap.png', dpi=300, bbox_inches='tight')


print("Received input files: ", sys.argv[1:])
germlineVariants=sys.argv[1]
somaticVariants=sys.argv[2]
subtractedVariants=sys.argv[3]
mutectVariants=sys.argv[4]
combinedCounts='/home/projects/dp_00005/data/nanbar/AnnotationScripts/Knijnenburg_Nasser_combined_countsOnly.txt'
makeHeatmap(germlineVariants, somaticVariants, subtractedVariants, mutectVariants, combinedCounts)

