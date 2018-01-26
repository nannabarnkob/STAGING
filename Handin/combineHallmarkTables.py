#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 15:39:52 2018

@author: nannabarnkob
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 15:26:00 2018

@author: nannabarnkob
"""
# resource for dataframe handling 
import pandas as pd
import re 
# for handling command-line input
import sys
import pdb

def combineDatasets(Knijnenburg, Nasser):
    # reading and preparing data 
    df_Knijnenburg = pd.read_table(Knijnenburg, index_col=0)
    df_Nasser = pd.read_table(Nasser, index_col=0)
    # print(df_Knijnenburg.index) 
    # print(df_Nasser.index)
    df_Knijnenburg = df_Knijnenburg.drop(['Other'])
    as_list = df_Nasser.index.tolist()
    idx = as_list.index('Resisting Cell Death ')
    as_list[idx] = 'Resisting Cell Death'
    df_Nasser.index = as_list
    hallmarks = df_Knijnenburg.index
    title = re.search('(.*)_HMmatches_Knijnenburg.txt',str(Knijnenburg))
    outfile = open(title.group(1)+'_HMmatches_combined.txt','w')
    outfile.writelines('Hallmark \t Count \t Gene symbols \n')
    kni_symbols_list = []
    nas_symbols_list = []
    all_kni = [] 
    all_nas = []
    for hallmark in hallmarks[0:len(hallmarks)-3]: 
        kni_symbols = df_Knijnenburg.loc[hallmark, ' Gene symbols ']
        #print(type(kni_symbols))
        if isinstance(kni_symbols, float):
            pass 
        else: 
            kni_symbols_list = kni_symbols.split(',')
        nas_symbols = df_Nasser.loc[hallmark, ' Gene symbols ']
        if isinstance(nas_symbols, float):
            pass 
        else: 
            nas_symbols_list = nas_symbols.split(',')
        all = nas_symbols_list + kni_symbols_list
        all_kni += kni_symbols_list
        all_nas += nas_symbols_list  
        #pdb.set_trace()
        outfile.writelines([hallmark, '\t', str(len(set(all))),'\t', ','.join(map(str, list(set(all)))),'\n'])
    sum_kni_nas = all_kni + all_nas
    outfile.writelines(['Total unique genes mapped to hallmarks \t',str(len(set(sum_kni_nas))),'\n'])
    outfile.writelines(['Number of gene mappings \t',str(len(sum_kni_nas)),'\n'])

Knijnenburg=sys.argv[1]
Nasser=sys.argv[2]
combineDatasets(Knijnenburg, Nasser)
