# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 23:12:27 2024

@author: imraa
"""

# Warning: All of the *tig*.gbk files, I gave them the prefix of the genome's name.
# So tig0000001.region003.gbk for WMMA2056 was manually renamed as file: WMMA2056.tig0000001.region003.gbk (example).

#%% Imports

# conda environment: python-class

# AntiSMASH relevant imports
import os
from bs4 import BeautifulSoup # beautiful soup version 4.10.0
from copy import deepcopy 
from os import chdir, getcwd
import re
import sys

import pandas as pd # pandas version 1.3.5
import numpy as np # numpy version 1.21.6
import matplotlib.pyplot as plt
import seaborn as sns

# BiG-SLiCE relevant imports
import sqlite3 # sqlite version 3.38.5
import pandas as pd # pandas version 1.3.5
import numpy as np # numpy version 1.21.6
# import sqlalchemy
import matplotlib.pyplot as plt # matplotlib version 3.5.1
# from matplotlib import rcParams
import math
import seaborn as sns # seaborn version 0.11.2
# import networkx as nx
import sklearn # scikit-learn version 1.0.2
# from sklearn import metrics
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
# import graphviz
# import pygraphviz
from sklearn.manifold import TSNE
from sklearn.manifold import MDS
import scipy # scipy version 1.7.3
import plotly.express as px # plotly express version 0.4.1
import plotly.io as pio # plotly version 5.10.0
# from sklearn.cluster import KMeans, AffinityPropagation, AgglomerativeClustering, Birch, DBSCAN, MeanShift, OPTICS # k-means clustering
# from sklearn.manifold import TSNE # t-SNE clustering
# from sklearn.model_selection import train_test_split # training/testing splitter
from scipy.cluster.hierarchy import linkage, dendrogram # Hierarchical Clustering
from numpy import unique, where
# from decimal import Decimal
# from pandas import ExcelWriter
# from sklearn import tree
# from mpl_toolkits import mplot3d
from statannotations.Annotator import Annotator
import copy
import plotly.express as px
import plotly.io as pio

# Dash
import dash
from dash import dcc, html, Input, Output, callback

import time
import pickle

#%% Existing Code

startTime = time.time();

antiSMASHdirectory = "C:/Users/imraa/Documents/UWM/BugniLab/Genome Assembly/antiSMASH/20240321-tim_request-genera";
scriptOutputs='??'

#%% Dictionary of antiSMASH product types converted to BiG-SCAPE classes

# Removed Siderophore as a category (they are now Others)    

bgSCAPE_conversion  = {
    "t1pks": "PKS I",
    "T1PKS": "PKS I",
    "transatpks": "PKS Other",
    "t2pks": "PKS Other",
    "t3pks": "PKS Other",
    "otherks": "PKS Other",
    "hglks": "PKS Other",
    "transAT-PKS": "PKS Other",
    "transAT-PKS-like": "PKS Other",
    "T2PKS": "PKS Other",
    "T3PKS": "PKS Other",
    "PKS-like": "PKS Other",
    "hglE-KS": "PKS Other",
    "prodigiosin": "PKS Other",
    "nrps": "NRPS",
    "NRPS": "NRPS",
    "NRPS-like": "NRPS",
    "thioamide-NRP": "NRPS",
    "NAPAA": "NRPS",
    "lantipeptide": "RiPP",
    "thiopeptide": "RiPP",
    "bacteriocin": "RiPP",
    "linaridin": "RiPP",
    "cyanobactin": "RiPP",
    "glycocin": "RiPP",
    "LAP": "RiPP",
    "lassopeptide": "RiPP",
    "sactipeptide": "RiPP",
    "bottromycin": "RiPP",
    "head_to_tail": "RiPP",
    "microcin": "RiPP",
    "microviridin": "RiPP",
    "proteusin": "RiPP",
    "guanidinotides": "RiPP",
    "lanthipeptide": "RiPP",
    "lipolanthine": "RiPP",
    "RaS-RiPP": "RiPP",
    "fungal-RiPP": "RiPP",
    "thioamitides": "RiPP",
    "lanthipeptide-class-i": "RiPP",
    "lanthipeptide-class-ii": "RiPP",
    "lanthipeptide-class-iii": "RiPP",
    "lanthipeptide-class-iv": "RiPP",
    "lanthipeptide-class-v": "RiPP",
    "ranthipeptide": "RiPP",
    "redox-cofactor": "RiPP",
    "RRE-containing": "RiPP",
    "epipeptide": "RiPP",
    "cyclic-lactone-autoinducer": "RiPP",
    "spliceotide": "RiPP",
    "crocagin": "RiPP",
    "RiPP-like": "RiPP", # was missing from BiG-SCAPE's class document (https://github.com/medema-group/BiG-SCAPE/wiki/big-scape-classes#hybrids)
    "amglyccycl": "Saccharide",
    "oligosaccharide": "Saccharide",
    "cf_saccharide": "Saccharide",
    "saccharide": "Saccharide",
    "terpene": "Terpene",
    "acyl_amino_acids": "Other",
    "arylpolyene": "Other",
    "aminocoumarin": "Other",
    "ectoine": "Other",
    "butyrolactone": "Other",
    "nucleoside": "Other",
    "melanin": "Other",
    "phosphoglycolipid": "Other",
    "phenazine": "Other",
    "phosphonate": "Other",
    "other": "Other",
    "cf_putative": "Other",
    "resorcinol": "Other",
    "indole": "Other",
    "ladderane": "Other",
    "PUFA": "Other",
    "furan": "Other",
    "hserlactone": "Other",
    "fused": "Other",
    "cf_fatty_acid": "Other",
    "siderophore": "Other", # Siderophore is its own category for now.
    "blactam": "Other",
    "fatty_acid": "Other",
    "PpyS-KS": "Other",
    "CDPS": "Other",
    "betalactone": "Other",
    "PBDE": "Other",
    "tropodithietic-acid": "Other",
    "NAGGN": "Other",
    "halogenated": "Other",
    "pyrrolidine": "Other",
    "mycosporine-like": "Other",
    "TfuA-related": "RiPP",
    # These later ones are necessary for antiSMASH v7
    'NI-siderophore': 'Other', # Was previously named Siderophore
    '2dos':'Saccharide', # aminoglycoside/aminocyclitol is considered a Saccharide, so I would consider 2-deoxy-streptamine aminoglycoside a Saccharide.
    'hydrogen-cyanide':'Other', # I don't think any other category fits.
    'NRP-metallophore':'NRPS', # Every time it appears, it's right on top of an annotated NRPS.
    'HR-T2PKS':'PKS Other', # Still a type 2 PKS.
    'aminopolycarboxylic-acid':'Other', # Labeled as other when it appears by itself in WMMD406.
    'phosphonate-like':'Other', # Assumed other.
}

#%% AntiSMASH analysis:

chdir(antiSMASHdirectory); # change the directory the antiSMASH directory specified previously.
wd = getcwd(); # get the current path, it should match the antiSMASh directory
# check if the directories match
if wd == antiSMASHdirectory:
    print('Currently in the antiSMASH directory...')

asmStorage = {};
htmlStorageFull = []; # container for the html file path names
htmlStorage = []; # container for mini html file path names
genomeNameS = []; # container for the individual genome names
titleStorage = [];

for filename in os.listdir(wd): # look at each folder/file within the main folder
    newF = os.path.join(wd, filename); # store the folder/file name, append, so we can access further
    genomeNameS.append(newF); # store genome names "generally speaking" here.
    for filename2 in os.listdir(filename): # look at each folder/file in the subfolder
        if filename2.endswith('x.html'): # find only the index.html file 
            fname = os.path.join(newF, filename2) # find the filename/path
            htmlStorageFull.append(fname); # store the file name path
            htmlStorage.append(os.path.join(filename, filename2));
            # print("Current file name ..", os.path.abspath(fname)) # print statement to see if we found the right files

# antiSMASH v7-specific processing (I think?)

# Scrape BGC names for later:
bgc_names = {}; # dictionary of dataframes
counter = 1;
for filename in os.listdir(wd): # look at each folder/file within the main folder
    newF = os.path.join(wd, filename); # store the folder/file name, append, so we can access further
    bgc_namesDF = pd.DataFrame();
    for filename2 in os.listdir(filename): # look at each folder/file in the subfolder
        if (('region' in filename2) and ('regions.js' not in filename2)): # find only the files containing 'region'
            bgc_namesDF = bgc_namesDF.append({'BGC Name':filename2.split('.gbk')[0]}, ignore_index = True)
    # bgc_names[transform_string(filename)] = bgc_namesDF;
    bgc_names[filename] = bgc_namesDF;
        
counter = 0; # useful counter
print('Iterating through ' + str(len(htmlStorage)) + ' folders in the antiSMASH directory...')
for d in htmlStorage: # for each html file found
    htmlfile = open(d) # don't forget to close it at the very end when we're done with it
    contents = htmlfile.read(); # read it out
    beautifulSoupT = BeautifulSoup(contents, 'html.parser') # can update to a later parser eventually
    titleStorage.append(beautifulSoupT.head.title); # this contains the # of regions found, might be useful
    
    regList = list(beautifulSoupT.find_all('tbody')); # only take the first one
    regList = list(regList[0]); # this contains all the info we could want
    
    asmStorageT = []; # storage for a singular html files import values (all regions)
    
    for i in np.arange(1, len(regList), 2): # 0, 2, 4,... all are empty headers '\n'
        # print(regList[i])
        poolT = regList[i].text.splitlines(); # get only the names, remove the \n (keeps spaces and +)
        asmStorageT.append(poolT); # take the split characters, put the list into the bigger list
    genomeInfo = pd.DataFrame(asmStorageT); # convert the big list into a dataFrame
    # this section is hardcoded, very bad idea, but works for now
    genomeInfo = genomeInfo.drop(columns = [0,1,3,4,6]); # remove the columns that are consistently bad under my current naming scheme
    genomeInfo.columns = ['Region', 'Type', 'From', 'To', 'Most Similar Known Cluster', 'MIBiG Type', 'Similarity'];
    nameP = genomeNameS[counter].split('\\')[-1]; # Assuming that the folder name is: '\path\to\folder\genomeName', this will take just the 'genomeName'
    # nameUpd = transform_string(nameP);
    nameUpd = nameP;
    asmStorage[nameUpd] = genomeInfo; # adds to dict the key (genome name), and the dataframe (genome info)
    htmlfile.close(); # close the file
    counter += 1
      
print('Finished iterating through the ' + str(counter) + ' folders in the antiSMASH directory, see asmStorage dictionary for output...')
print('Final amount of antiSMASH results used: ' + str(len(asmStorage)))
asmStorageTemp = asmStorage.copy();

#%% Convert BGC product type to BiG-SCAPE classes:

print('Replacing antiSMASH product types with BiG-SCAPE classes...')
asmStorageClean = deepcopy(asmStorageTemp);
keyList = list(asmStorageClean.keys());
counter1 = 0;
# Replace antiSMASH product types with BiG-SCAPE classifications.
    # asmStorageClean will contain the replaced information (replaced inplace)
    # asmStorageTemp contains the old antiSMASH product types
for i in keyList: # access each dataframe
    # print(i);
    counter1 +=1;
    # print(counter1);
    dFClean = asmStorageClean[i];
    dfClean_fix = dFClean['Type'].apply(lambda x: ', '.join([bgSCAPE_conversion[word] for word in x.split(',')]));
    dFClean.Type = dfClean_fix;

# Get all relevant antiSMASH info, stored here. (Unnecessary?)
typeDF_List = [];
for key, df in asmStorageClean.items():
    # print(key)
    df['Genome'] = key;
    if key != 'WMMC500': # Hardcoded because I don't want WMMC500
        typeDF_List.append(df)
combinedTypeDF_List = pd.concat(typeDF_List, ignore_index=True)

# Map the BGC names to the Product Types:
bgcClassStorage = {};
for key, df in asmStorageClean.items():
    mixStorageDF = pd.DataFrame({'name 1': bgc_names[key]['BGC Name'].values, 'bgc product': df.Type.values.tolist()});
    bgcClassStorage[key] = mixStorageDF;

bgSCAPEclasses = pd.concat(bgcClassStorage, ignore_index = True);

#%% Scheme for devising Hybrid (not PKS/NRPS) category:

def clean_and_merge(s):
    return ', '.join(sorted(set(map(str.strip, s.split(',')))))

# Rule: If there's an instance of 'X, X, Y', consider it 'X, Y'. This applies for any number of X, and so on.
bgSCAPEclasses_copy = deepcopy(bgSCAPEclasses);
bgSCAPEclasses_copy['bgc product'] = bgSCAPEclasses_copy['bgc product'].apply(clean_and_merge);

def check_hybrid_notPKS_NRPS(bgc_product):
    # If the specific condition for PKS/NRPS is found.
    if bgc_product == 'PKS I, NRPS' or bgc_product == 'NRPS, PKS I':
        return 'PKS/NRPS'
    # If there are 2 or more elements in the bgc product (separated by a string), and they're not the PKS/NRPS.
    elif ',' in bgc_product:
        return 'Hybrid (Not PKS/NRPS)'
    # If the previous conditions have failed because there is no PKS/NRPS and there's only one element.
    return bgc_product

# Rule: If not 'PKS I, NRPS' or vice-versa, and not one element, then label with 'Hybrid (Not PKS/NRPS)'.
bgSCAPEclasses_copy2 = deepcopy(bgSCAPEclasses_copy);
bgSCAPEclasses_copy2['bgc product'] = bgSCAPEclasses_copy2['bgc product'].apply(check_hybrid_notPKS_NRPS);

# Calculations:

# Extract genome name from 'name 1' column
bgSCAPEclasses_copy2['genome_name'] = bgSCAPEclasses_copy2['name 1'].apply(lambda x: x.split('.')[0])
# Group by 'genome_name' and 'bgc product', then count occurrences
bgc_counts = bgSCAPEclasses_copy2.groupby(['genome_name', 'bgc product']).size().reset_index(name='count')
# Pivot the table to get 'bgc product' as columns
bgc_counts_pivot = bgc_counts.pivot(index='genome_name', columns='bgc product', values='count').fillna(0)
# Reset index to make 'genome_name' a column instead of index
bgc_counts_pivot.reset_index(inplace=True)
# Reorder the columns
desired_order = ['genome_name','Terpene', 'PKS I', 'RiPP', 'PKS Other', 'PKS/NRPS', 'NRPS', 'Other', 'Saccharide', 'Hybrid (Not PKS/NRPS)']
bgc_counts_pivot_final = bgc_counts_pivot.reindex(columns=desired_order)
# Sum across rows, ignoring the 'genome_name' column
bgc_counts_pivot_final['Total'] = bgc_counts_pivot_final.drop(columns=['genome_name']).sum(axis=1)

# Can check bgc_counts_pivot_final to see the BGCs.




