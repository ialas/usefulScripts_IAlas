# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 21:18:42 2023

@author: imraan alas
"""

# The webpage will be at localhost:8002, or 127.0.0.1:8002

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

# AntiSMASH folders/files:
antiSMASHdirectory = input('What is the filepath to the folder containing the AntiSMASH data? If nothing is typed, will assume default.\nDirectory: ')
if os.path.exists(antiSMASHdirectory):
    print('AntiSMASH version 5.1.1 directory: ' + antiSMASHdirectory)
else:
    antiSMASHdirectory = r"C:/Users/imraa/Documents/UWM/BugniLab/Genome Assembly/micromonosporaceae-AntiSmash511-Clean";
    print('Could not find antiSMASH directory, using default directory:\n' + antiSMASHdirectory)

# BiG-SLiCE folders/files:
pathToDataDB = input('Filepath to the data.db for the already run BiG-SLiCE data?\nFilepath: ')
if os.path.exists(pathToDataDB):
    print('The path to the data.db file in the result folder: ' + pathToDataDB)
else:
    pathToDataDB = "C:/Users/imraa/Downloads/DataS10-bigsliceOutputs-antiSMASH-v511/DataS10-bigsliceOutputs-antiSMASH-v511/full_run_result/result/data.db"; # This is the path to the data.db file, once copied from Linux and put into Windows.
    print('Couldn\'t find the data.db file, here\'s an example location: ' + pathToDataDB)
    
pathToReportsDB = input('Filepath to the reports folder for the already run BiG-SLiCE data?\nFilepath: ')
if os.path.exists(pathToReportsDB):
    print('The path to the reports folder:\n' + pathToReportsDB)
else:
    pathToReportsDB = "C:/Users/imraa/Downloads/DataS10-bigsliceOutputs-antiSMASH-v511/DataS10-bigsliceOutputs-antiSMASH-v511/full_run_result/result/data.db"; # This is the path to the data.db file, once copied from Linux and put into Windows.
    print('Couldn\'t find the reports folder, here\s an example location: ' + pathToReportsDB)

run_id = int(input('Choose a run_id that was used for querying your BGCs against BiG-SLiCE\'s pre-processed database.\n1 is using no Threshold.\n2 is using Threshold 300.\n'+
               '4 is using Threshold 600.\n6 is using Threshold 900.\n7 is using Threshold 1200.\n8 is using Threshold 1500.\nEnter run_id here: '))
if run_id not in [1,2,4,6,7,8]:
    print('ERROR: Please select a proper run_id, with numbers possible either 1, 2, 4, 6, 7, 8')
    sys.exit()
else:
    runThresholdPair = {1: 0, 2: 300, 4: 600, 6: 900, 7: 1200, 8: 1500};
    thresholdValue = runThresholdPair[run_id];
    print('You have selected run_id ' + str(run_id) + ' with associated threshold value of ' + str(thresholdValue))

reportsRuns = input('Filepath to a text file where each line contains a number that represents the BiG-SLiCE runs used for this analysis.\nFor example, if you wanted to run reports 100-102, the text file would contain:\n100\n101\n102\nFilepath: ')
if os.path.exists(reportsRuns):
    # Read the folder list from the file
    with open('reportsRun', 'r') as file:
        # Assuming each line contains a folder number
        allFolders = [int(line.strip()) for line in file]
        print('The report folders used for this analysis are: ')
        print(allFolders)
else:
    print('ERROR: Couldn\'t find file containing the BiG-SLiCE runs intended to be analyzed...')
    print('Reminder, this should be a simple text file where each line is a number that represents the run you wish to include in the analysis')
    # Default case (for my own files):
    allFolders = np.concatenate([np.arange(111,145,1), np.arange(149,157,1)]); # antismash v5.1.1: Right now, it looks at reports folders 111 to 144, and 150 to 156. The numbers can be found in the BiG-SLiCE visualization.
    print('(DEFAULT): The report folders used for this analysis are: ')
    print(allFolders)

# A little bit of pre-processing for my own data.
if antiSMASHdirectory == r"C:/Users/imraa/Documents/UWM/BugniLab/Genome Assembly/micromonosporaceae-AntiSmash511-Clean":
    strainsToRemoveB = "a1363|b482|b486"; # I remove these strains in the BiG-SLICE analysis because they failed QC.
    print('Strains to remove from the BiG-SLiCE visualizations: ' + strainsToRemoveB)

# Determine script output folder:
scriptOutputs = input('What directory should files generated by this analysis be stored?\nDirectory:')
if os.path.exists(scriptOutputs):
    print('Relevant outputs will be stored in: ' + scriptOutputs)
else:
    scriptOutputs = 'C:/Users/imraa/Documents/UWM/BugniLab/Genome Assembly/scriptOutputs'
    print('The directory does not appear to exist. Here\'s an example one:\n' + scriptOutputs)

# Check for already run files.
print('If you\'ve run this analysis already, this section double-checks to see if there\'s any stored files it can use.')
# Path to csv file containing BGCs x GCFs (BiG-SLiCE): Used to minimize resource usage after running Section 5 once.
if os.path.isfile(scriptOutputs + '/' + 'antismash5-strainsQueried-rank.csv'):
    queriedBGC_GCFpath = scriptOutputs + '/' + 'antismash5-strainsQueried-rank.csv';
    print('Found csv file containing BGCs x GCFs, here:\n' + queriedBGC_GCFpath)
else:
    print('No previously analyzed data regarding BiG-SLiCE queries run through this program already exists, will run and store in:\n' + queriedBGC_GCFpath)

if os.path.isfile(scriptOutputs + '/' + 'mds_data_run.pickle'):
    mdsFile = scriptOutputs + '/' + 'mds_data_run.pickle';
    print('Found pickled file containing MDS results, here:\n' + mdsFile)
else:
    print('Could not find already run MDS data, will generate independently and save into ' + scriptOutputs + '/mds_run_data.pickle')

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

# antiSMASH v5.1.1 specific processing

# Fix the "a998" to "WMMA998":
def transform_string(input_str):
    # Extract the first letter
    first_letter = input_str[0].upper()
    # Add the relevant stuff & replace the first letter.
    result_str = "WMM" + first_letter + input_str[1:]
    return result_str
# Scrape BGC names for later:
bgc_names = {}; # dictionary of dataframes
counter = 1;
for filename in os.listdir(wd): # look at each folder/file within the main folder
    newF = os.path.join(wd, filename); # store the folder/file name, append, so we can access further
    bgc_namesDF = pd.DataFrame();
    for filename2 in os.listdir(filename): # look at each folder/file in the subfolder
        if (('region' in filename2) and ('regions.js' not in filename2)): # find only the files containing 'region'
            bgc_namesDF = bgc_namesDF.append({'BGC Name':filename2.split('.gbk')[0]}, ignore_index = True)
    bgc_names[transform_string(filename)] = bgc_namesDF;
        
counter = 0; # useful counter
print('Iterating through ' + str(len(htmlStorage)) + ' folders in the antiSMASH directory...')
for d in htmlStorage: # for each html file found
    htmlfile = open(d) # don't forget to close it at the very end when we're done with it
    contents = htmlfile.read(); # read it out
    beautifulSoupT = BeautifulSoup(contents, 'html.parser') # can update to a later parser eventually
    titleStorage.append(beautifulSoupT.head.title); # this contains the # of regions found, might be useful
    
    regList = list(beautifulSoupT.find_all('tbody')); # only take the first one
    regList = list(regList[-1]); # this contains all the info we could want
    
    asmStorageT = []; # storage for a singular html files import values (all regions)
    
    for i in np.arange(1, len(regList), 2): # 0, 2, 4,... all are empty headers '\n'
        poolT = regList[i].text.splitlines(); # get only the names, remove the \n (keeps spaces and +)
        asmStorageT.append(poolT); # take the split characters, put the list into the bigger list
    genomeInfo = pd.DataFrame(asmStorageT); # convert the big list into a dataFrame
    # this section is hardcoded, very bad idea, but works for now
    genomeInfo = genomeInfo.drop(columns = [0,1,3,4,6,9]); # remove the columns that are consistently bad under my current naming scheme
    genomeInfo.columns = ['Region', 'Type', 'From', 'To'];
    nameP = genomeNameS[counter].split('\\')[-1]; # Assuming that the folder name is: '\path\to\folder\genomeName', this will take just the 'genomeName'
    nameUpd = transform_string(nameP);
    asmStorage[nameUpd] = genomeInfo; # adds to dict the key (genome name), and the dataframe (genome info)
    htmlfile.close(); # close the file
    counter += 1
    
# If there was any contamination that you are able to identify in your genome, and want to filter it out from the data here.
    # Example: If there's contamination in my bacteria strain d975, and I found that it was in antiSMASH regions 25 through 38 (counting down from the top in antismash), this is how I would remove those data pieces.
print('The code will now attempt to remove extraneous elements in my data...\nIf they are not present, you will see some print statements indicating that...')
# Remove stuff from asmStorage
try:
    asmStorage['WMMD975'].drop(list(np.arange(25,39,1)),inplace=True); # Because python is 0 indexed, region 1 is actually the 0th region in python.
except:
    print('WMMD975 extraneous information was not present...')
try:
    asmStorage.pop("WMMA1363");
except:
    print('WMMA1363 was not present...')
try:
    asmStorage.pop("WMMB482");
except:
    print('WMMB482 was not present...')
try:
    asmStorage.pop("WMMB486");
except:
    print('WMMB486 is not present...')
try:
    asmStorage.pop('WMMC500');
except:
    print('WMMC500 is not present...')
# Remove stuff from bgc_names
try:
    bgc_names.pop("WMMA1363")
except:
    pass
try:
    bgc_names.pop("WMMB482")
except:
    pass
try:
    bgc_names.pop("WMMB486")
except:
    pass
try:
    bgc_names.pop('WMMC500')
except:
    pass
try:
    bgc_names['WMMD975'].drop(list(np.arange(25,39,1)),inplace=True); # Because python is 0 indexed, region 1 is actually the 0th region in python.
except:
    pass
        
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

bgSCAPEclasses = pd.concat(bgcClassStorage, ignore_index = True); # matches the old bgSCAPEclasses excel document

    
#%% BiG-SLiCE (Data Generation): Unnecessary to run in future attempts.

try:
    chdir(scriptOutputs)
    # See if the MDS data is already present
    # Load pickle: (contains already analyzed MDS data)
    with open(mdsFile, 'rb') as handle:
        dictDFMDS_ALL = pickle.load(handle)
except:
    chdir(scriptOutputs)
    # Check to see if the work has been done already.
    try:
        superDF3 = pd.read_csv(queriedBGC_GCFpath);
        print('Taking already run BiG-SLiCE analysis...')
    except:
        print('Running full BiG-SLiCE analysis...')
    
        # Computationally-intensive. (Generates a BGC x GCF dataframe)
        connMib = sqlite3.connect(pathToDataDB);
        currMib = connMib.cursor();
        # allFolders = np.concatenate([np.arange(111,145,1), np.arange(149,157,1)]); # antismash v5.1.1 # DEFINED IN SECTION 3
        storageDict = {};
        counterN=0;
        for j in allFolders: # whatever the report runs i want to look at are
            print('Analyzing folder: ' + str(j) + ' and scraping all queried distances associated with this strain...')
            connF = sqlite3.connect(pathToReportsDB + str(j) + "/data.db") # open a connection
            curF = connF.cursor();
            dfGCF = pd.read_sql_query("SELECT * FROM gcf_membership", connF);
            dfSliceGCF = dfGCF[(dfGCF.iloc[:,3] == 0)]; # get only the highest rank values
            listBGC = list(dfSliceGCF.gcf_id);
            numBGCs = len(dfGCF.bgc_id.value_counts());
            list1, list2, list3, list4 = list(), list(), list(), list();
            count12 = 0;
            for i in np.arange(1,numBGCs+1,1):
                name1 = curF.execute("SELECT name FROM bgc WHERE id=" + str(i)).fetchall()[0]; # modified for this new version
                dist1 = curF.execute("SELECT membership_value, gcf_membership.gcf_id FROM gcf_membership WHERE bgc_id=" + str(i)).fetchall(); # modded to grab all ranks
                # gcf1 = currMib.execute("SELECT id_in_run FROM gcf, clustering WHERE gcf.clustering_id=4 AND gcf.id=" + str(dist1[1]) + " AND clustering.run_id=4").fetchall()[0];
                name2 = name1[0].split('/')[0]; # genome name folder
                name22 = name1[0].split('/')[1]; # tig.region
                # I HARDCODED THIS SECTION because it was easier. 
                # What this does is it scrapes the genome name associated with this run, and then the contig and the region for the BGC. You can re-implement this how you see fit.
                if name2 != 'm6_a1363': # hard code it because i messed up.
                    name3 = name2.split('.')[-2]; # bc_genomeName
                    if len(name3.split('_')) < 2: # a rule for if its bc12_genomeName or bc-genomeName
                        name4 = name3.split('-')[-1];
                        name42 = name4 + '.' + name22;
                    else:
                        name4 = name3.split('_')[-1]; # just genomeName
                        name42 = name4 + '.' + name22; # genomeName.tig00.region00
                else:
                    name4 = 'a1363';
                    name42 = name4 + '.' + name22;
                for k in np.arange(0,len(dist1),1):
                    gcf1 = currMib.execute("SELECT id_in_run FROM gcf, clustering WHERE gcf.clustering_id=4 AND gcf.id=" + str(dist1[k][1]) + " AND clustering.run_id=4").fetchall()[0];
                    list1.append(name42);
                    list2.append(dist1[k][0]);
                    list3.append(dist1[k][1]);
                    list4.append(gcf1[0]);
                count12 += 1;
            counterN +=1;
            print('Completed folder: ' + str(counterN) + ' of ' + str(len(allFolders)))
            genomeName = name4; # name1[0].split('.')[0];
            nameDistDict = {'BGC': list1, 'Distance': list2, 'gcf_ID': list3, 'GCF_Value':list4};
            nameDist = pd.DataFrame(nameDistDict);
            storageDict[genomeName] = nameDist;
            connF.close();
        
        # REMOVE IF ERROR (Hardcoded for my own data)
        # Remove the extra stuff in WMMD975 that might be contamination.
        try:
            d975Stuff = storageDict['d975'];
            d975StuffKeep = d975Stuff[d975Stuff['BGC'].str.contains('tig00000001')];
            d975StuffRemove = d975Stuff[~d975Stuff['BGC'].str.contains('tig00000001')];
            storageDict['d975'] = d975StuffKeep;
        except:
            print('WMMD975 was not found in the BiG-SLiCE files...')
        
        print('Finished scraping all of the BiG-SLiCE report folders...')
        print('Construcing file for exporting...')
        # Merge all of the information into one big dataframe. (X, where X = GCF*BGC, by 4)
        newDF = pd.DataFrame(columns = list(storageDict['a1363'].columns));
        dfs_list = [];
        for i in list(storageDict.keys()):
            print(i)  # Tracker: outputs to terminal
            oldDF = storageDict[i].sort_values(by=['BGC', 'gcf_ID'])
            dfs_list.append(oldDF)
        newDF = pd.concat(dfs_list, ignore_index=True)
        testDF = newDF;
        
        # Dataframe (rows are BGCs, columns are GCF values)
        testT1 = list(testDF.BGC); # get the BGC names.
        testT1 = list(set(testT1)); # get the unique BGC names.
        testT2 = list(testDF.gcf_ID); # get the gcf id
        testT2 = list(set(testT2)); # get the unique GCF ids
        
        # Get the individual data into lists, ordered and sorted.
        superList1 = [];
        counterT = 0;
        counterT2 = 0;
        for i in testT1:
            placeHolder1 = testDF.loc[(testDF['BGC'] == i)]; # returns rows of only the BGC (sorted by gcf_ID)
            pH2 = list(placeHolder1.Distance);
            pH3 = list(placeHolder1.gcf_ID);
            superList1.append(pH2)
            counterT2 += 1;
            if (counterT2 - counterT) >= 30: # Tracker (outputs to console)
                # print(counterT2);
                counterT = counterT2;
        
        print("Making a mega dataframe, takes awhile...")
        # Make a super dataframe, with the columns as GCFs-200946 to the end.
        superDF2 = pd.DataFrame([i for i in superList1], columns = testT2);
        superDF2.index = testT1; # change the row names to be the BGC names.
        
        superDFbackup = superDF2.copy(deep=True);
        
        # REMOVE IF ERROR, BASED ON MY ORIGINAL DATA: Removes the Streptomycetaceae (c500) and the ones that failed QC.
        superDF2 = superDF2[~superDF2.index.str.contains('c500|a1363|b482|b486')];
        
        # print('Exporting the BGC x GCF file for future re-use to ' + resultsDirectory + r'\antismash5-strainsQueried-rank.csv')
        
        # Must swap the rows and columns, as excel doesn't like having 29,955 columns.
        superDF2T = superDF2.T;
        
        # Export the superDF2 (GCF x BGC file) to CSV for ease of access later.
        if not os.path.isfile('antismash5-strainsQueried-rank.csv'):
            superDF2T.to_csv('antismash5-strainsQueried-rank.csv');
        
    #%% BiG-SLiCE (can run without previous section, if it's been run before):
    
    chdir(scriptOutputs)    
    
    # ---- If re-using a csv file containing all BGCs x GCFs queried, can load here.
    print('Loading file at ' + queriedBGC_GCFpath)
    # Load file.
    try:
        superDF3 = pd.read_csv(queriedBGC_GCFpath);
    except:
        print('You need to run the previous section, titled:\nBiG-SLiCE (Data Generation)')
    # Get the index.
    superDF3.index = superDF3.iloc[:,0]; # replace index with GCFs
    # Remove the first column (which previously contained the index)
    superDF3.drop(columns = superDF3.columns[0], axis=1, inplace=True); # remove first column
    # Transpose so it is in the shaep of BGC x GCF.
    superDF3=superDF3.T;
    # Scrape the index.
    testT1 = list(superDF3.index);
    
    superDF2 = superDF3;
    
    # ---- Data necessary for the Figures.
    
    #%% Generate Data for the non-Hybrid version
    
    # Load bigscape class product types in relation to bgc names.
    # bgSCAPEclasses = pd.read_excel(bigSCAPEclassesAllpath);
    bgSCAPEclasses = bgSCAPEclasses
    
    # Create filters for eventually splitting into the relevant datasets:    
    
    # List of substrings to check
    substrings = ['Terpene', 'PKS I', 'RiPP', 'Other', 'PKS Other', 'PKS/NRPS', 'NRPS', 'Saccharide']
    # Create a new column for each category and mark with True or False
    for sub in substrings:
        bgSCAPEclasses[sub] = bgSCAPEclasses['bgc product'].str.contains(sub, case=False, na=False, regex=True)
    # Create a new column for 'PKS/NRPS'
    bgSCAPEclasses['PKS/NRPS'] = bgSCAPEclasses['bgc product'].str.contains('PKS I') & bgSCAPEclasses['bgc product'].str.contains('NRPS')
    bgSCAPEclasses['PKS Other'] = bgSCAPEclasses['bgc product'].str.contains(r'PKS Other\b', case=False, na=False, regex=True)
    bgSCAPEclasses['Other'] = bgSCAPEclasses['bgc product'].str.contains(r'(?<!PKS )Other', case=False, na=False, regex=True)
    
    # Filter the DataFrame based on the presence of substrings
    filtered_df = bgSCAPEclasses[bgSCAPEclasses[substrings].any(axis=1)]
    missing_df = bgSCAPEclasses[~bgSCAPEclasses[substrings].any(axis=1)]
    
    colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#e31a1c','#fdbf6f','#ff7f00','#cab2d6']; # colors match the ITOL TREE
    classes1 = ['Terpene', 'PKS I', 'RiPP', 'PKS Other', 'PKS/NRPS', 'NRPS', 'Other', 'Saccharide']; # matches the colors indices.
    
    # Construct a dictionary of palettes.
    palDict = {};
    for i in range(0,len(colors),1):
        palDict[classes1[i]] = colors[i];
        
    classDict_Full = {};
    for i in classes1:
        classDict_Full[i] = bgSCAPEclasses[bgSCAPEclasses[i]];
    
    # # Find which indices in the BGCs map to which classes.
    # classDict = {}; # Keys are the classes, each key links to a list of the indices where it was found in classList
    # for i in classes1:
    #     classDict[i] = [a for a, x in enumerate(classList) if x == i]; 
    # Copy superDF2 just incase.
    superCatDF2 = superDF2.copy(deep=True);
    
    dictDistMetrics = {}; # initialize empty dictionary for storage
    distMets = ['euclidean', 'cosine', 'cityblock', 'l1', 'l2', 'manhattan', 'braycurtis', 'canberra', 'chebyshev', 'correlation', 'hamming'] # this is where the different strings for distance metrics go
    superDict0 = {}; # dictionary of the different simplified datasets I want to play with.
    # Create the datasets (sliced by bigSCAPE categories)
    for b in classes1: 
        # Store the relevant data from superDF2 based on the bigscape class the BGC contains.
        superDict0[b+"_orig"] = superCatDF2[superCatDF2.index.isin(classDict_Full[b]['name 1'])]
        
    # Initialize an empty dictionary
    dictDistMetrics = {};
    # Run the distance metrics, the structure of dictDistMetrics is as follows:
        # 'Terpene_orig' -> 'Terpene_orig_euclidean", "Terpene_orig_cosine", etc (each is a dataframe)
    for j in list(superDict0.keys()):
        dictDistMetrics[j] = {};
        print('Iterating over BGC class: ' + j)
        for i in distMets:
            print('Computing distance for: ' + i)
            outputDist0 = pd.DataFrame(sklearn.metrics.pairwise_distances(superDict0[j], Y=None, metric=i, n_jobs=None, force_all_finite=True), columns = superDict0[j].index.tolist());
            outputDist0.index = superDict0[j].index.tolist();
            dictDistMetrics[j][j+"_"+i] = outputDist0; # pop it back in
    # Convert the distance metrics from dataframes to arrays for MDS.
    print('Converting dataframes to arrays...')
    dictArray = {};          
    for j in list(superDict0.keys()):
        dictArray[j] = {};
        for i in distMets:
            outputArray = dictDistMetrics[j][j+"_"+i].to_numpy();
            dictArray[j][j+"_"+i] = outputArray;  
    
    # Perform MDS
    print('Performing MDS...')
    dictMDS = {};        
    nComps = 2; # number of components
    for j in list(superDict0.keys()):
        print('Iterating over BGC class: ' + j)
        dictMDS[j] = {};
        for i in distMets:
            print('Metric used for MDS currently: ' + i)
            mdsModelA = MDS(n_components = nComps, dissimilarity='precomputed', random_state=0);
            dataTrans = mdsModelA.fit_transform(dictArray[j][j+"_"+i]);
            dictMDS[j][j+"_"+i] = dataTrans;
    
    print('Converting to dataframes for ease of usage...')
    # aveDist = 3; # No longer necessary. I set hardcoded values of 1,2,3,4,5 as k in the below loop.
    dictDFMDS = {};
    for j in list(superDict0.keys()):
        dictDFMDS[j] = {};
        superDict0_j = superDict0[j];
        print('Iterating through: ' + j)
        for i in distMets:
            dfMDStemp = pd.DataFrame(dictMDS[j][j+"_"+i], index = superDict0_j.index.tolist())
            dfMDStemp['Category'] = j.split("_")[0]; # bigscape product types
            # grab the <900 >900 information from bigslice
            # list comprehension (if/else statement). GOal: '> 900' if the BGC fell above 900 in BiGSLICE, otherwise '< 900'
            # bigSLICEinfo = ['> 900' if storage2DF.loc[storage2DF['BGC'] == iter1].Distance.tolist()[0] >= 900 else '< 900' for iter1 in superDict0[j].index.tolist()];
            
            # Removed dependency on storageDict
            bigSLICEinfo_v2 = ['> 900' if min(superDF2.loc[iter1]) >= 900 else '< 900' for iter1 in superDict0_j.index.tolist()]
            # find the distance to the closest GCF
            # closestGCFDist = [storage2DF.loc[storage2DF['BGC'] == iter1].Distance.tolist()[0] for iter1 in superDict0[j].index.tolist()];
            
            # Removed dependency on storageDict
            closestGCFDist_v2 = [min(superDF2.loc[iter1]) for iter1 in superDict0_j.index.tolist()]
            
            # Do some average distance calculations.
            for k in [1,2,3,4,5]:
                # find the distance to the X closest GCFs, then average them. Where X is set above as aveDist
                # aveDistStorage = [sorted(superDict0_j.loc[iter2])[:k] for iter2 in superDict0_j.index]; # list of lists
                # averageDist = [np.mean(entry1) for entry1 in aveDistStorage]; # average each list inside the list of lists
                # dfMDStemp['Average of ' + str(k) + ' Distances'] = averageDist;
                dfMDStemp[f'Average of {k} Distances'] = [np.mean(sorted(superDict0_j.loc[iter2])[:k]) for iter2 in superDict0_j.index]
    
            
            dfMDStemp['BIG-SLICE'] = bigSLICEinfo_v2; # not dependent on storageDict anymore
            dfMDStemp['Distance'] = closestGCFDist_v2; # not dependent on storageDict anymore
            # dfMDStemp['Average Distance'] = averageDist;
            # dfMDStemp['# of GCF Distances used'] = aveDist;
            # Add metadata
            # bgSCAPEinfo = [bgSCAPEclasses.loc[bgSCAPEclasses['name 1'] == i]['bgc product'].values for i in dfMDStemp.index];
            # dfMDStemp['Metadata'] = bgSCAPEinfo
            dfMDStemp['Metadata'] = dfMDStemp.index.map(lambda x: bgSCAPEclasses.loc[bgSCAPEclasses['name 1'] == x]['bgc product'].values[0])
            dfMDStemp.rename(columns={0:'x', 1:'y'}, inplace=True) # clean up column names for ease of reference later
            
            dictDFMDS[j][j+"_"+i] = dfMDStemp;
    
    specificMets = distMets # this is where the different strings for distance metrics go
    dictDFMDS_ALL = {};
    dictDFMDS_ALL['dictDFMDS_noHybridCategory'] = dictDFMDS;
    
    #%% Generate data necessary for the Hybrid dataset.
    
    bgSCAPEclasses = pd.concat(bgcClassStorage, ignore_index = True); # matches the old bgSCAPEclasses excel document
    
    # Create filters for eventually splitting into the relevant datasets:    
    
    # List of substrings to check
    substrings = ['Terpene', 'PKS I', 'RiPP', 'Other', 'PKS Other', 'PKS/NRPS', 'NRPS', 'Saccharide']
    # Create a new column for each category and mark with True or False
    for sub in substrings:
        bgSCAPEclasses[sub] = bgSCAPEclasses['bgc product'] == sub;
        # bgSCAPEclasses['bgc product'].str.match(sub, case=True, na=False)
    # Create a new column for 'PKS/NRPS'
    pks_nrps_substrings = ['PKS I, NRPS', 'NRPS, PKS I'];
    bgSCAPEclasses['PKS/NRPS'] = [value in pks_nrps_substrings for value in bgSCAPEclasses['bgc product']]
    
    new_substrings = ['Terpene', 'PKS I', 'RiPP', 'Other', 'PKS Other', 'PKS/NRPS', 'NRPS', 'Saccharide']
    
    # Filter the DataFrame based on the presence of substrings
    filtered_df = bgSCAPEclasses[bgSCAPEclasses[new_substrings].any(axis=1)]
    missing_df = bgSCAPEclasses[~bgSCAPEclasses[new_substrings].any(axis=1)]
    
    # Currently, the Hybrid rule is this:
        # The BGCs that are in each category, they are purely just 1 BGC. Meaning, 'PKS I' only shows up in the category 'PKS I'.
        # So, 'PKS I, PKS I' would show up only in the Hybrid category. 'PKS I, Other' would show up only in the Hybrid category.
        # However, 'PKS/NRPS' category includes only BGCs that have 'PKS I, NRPS' or 'NRPS, PKS I'.
    # Create the Hybrid category.
    bgSCAPEclasses['Hybrid (Not PKS/NRPS)'] = ~bgSCAPEclasses[new_substrings].any(axis=1);
    
    colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a']; # colors match the ITOL TREE
    classes1 = ['Terpene', 'PKS I', 'RiPP', 'PKS Other', 'PKS/NRPS', 'NRPS', 'Other', 'Saccharide', 'Hybrid (Not PKS/NRPS)']; # matches the colors indices.
    
    # Construct a dictionary of palettes.
    palDict = {};
    for i in range(0,len(colors),1):
        palDict[classes1[i]] = colors[i];
        
    classDict_Full = {};
    for i in classes1:
        classDict_Full[i] = bgSCAPEclasses[bgSCAPEclasses[i]];
    
    # # Find which indices in the BGCs map to which classes.
    # classDict = {}; # Keys are the classes, each key links to a list of the indices where it was found in classList
    # for i in classes1:
    #     classDict[i] = [a for a, x in enumerate(classList) if x == i]; 
    # Copy superDF2 just incase.
    superCatDF2 = superDF2.copy(deep=True);
    
    dictDistMetrics = {}; # initialize empty dictionary for storage
    distMets = ['euclidean', 'cosine', 'cityblock', 'l1', 'l2', 'manhattan', 'braycurtis', 'canberra', 'chebyshev', 'correlation', 'hamming'] # this is where the different strings for distance metrics go
    superDict0 = {}; # dictionary of the different simplified datasets I want to play with.
    # Create the datasets (sliced by bigSCAPE categories)
    for b in classes1: 
        # Store the relevant data from superDF2 based on the bigscape class the BGC contains.
        superDict0[b+"_orig"] = superCatDF2[superCatDF2.index.isin(classDict_Full[b]['name 1'])]
        
    # Initialize an empty dictionary
    dictDistMetrics = {};
    # Run the distance metrics, the structure of dictDistMetrics is as follows:
        # 'Terpene_orig' -> 'Terpene_orig_euclidean", "Terpene_orig_cosine", etc (each is a dataframe)
    for j in list(superDict0.keys()):
        dictDistMetrics[j] = {};
        print('Iterating over BGC class: ' + j)
        for i in distMets:
            print('Computing distance for: ' + i)
            outputDist0 = pd.DataFrame(sklearn.metrics.pairwise_distances(superDict0[j], Y=None, metric=i, n_jobs=None, force_all_finite=True), columns = superDict0[j].index.tolist());
            outputDist0.index = superDict0[j].index.tolist();
            dictDistMetrics[j][j+"_"+i] = outputDist0; # pop it back in
    # Convert the distance metrics from dataframes to arrays for MDS.
    print('Converting dataframes to arrays...')
    dictArray = {};          
    for j in list(superDict0.keys()):
        dictArray[j] = {};
        for i in distMets:
            outputArray = dictDistMetrics[j][j+"_"+i].to_numpy();
            dictArray[j][j+"_"+i] = outputArray;  
    
    # Perform MDS
    print('Performing MDS...')
    dictMDS = {};        
    nComps = 2; # number of components
    for j in list(superDict0.keys()):
        print('Iterating over BGC class: ' + j)
        dictMDS[j] = {};
        for i in distMets:
            print('Metric used for MDS currently: ' + i)
            mdsModelA = MDS(n_components = nComps, dissimilarity='precomputed', random_state=0);
            dataTrans = mdsModelA.fit_transform(dictArray[j][j+"_"+i]);
            dictMDS[j][j+"_"+i] = dataTrans;
    
    print('Converting to dataframes for ease of usage...')
    # aveDist = 3; # No longer necessary. I set hardcoded values of 1,2,3,4,5 as k in the below loop.
    dictDFMDS = {};
    for j in list(superDict0.keys()):
        dictDFMDS[j] = {};
        superDict0_j = superDict0[j];
        print('Iterating through: ' + j)
        for i in distMets:
            dfMDStemp = pd.DataFrame(dictMDS[j][j+"_"+i], index = superDict0_j.index.tolist())
            dfMDStemp['Category'] = j.split("_")[0]; # bigscape product types
            # grab the <900 >900 information from bigslice
            # list comprehension (if/else statement). GOal: '> 900' if the BGC fell above 900 in BiGSLICE, otherwise '< 900'
            # bigSLICEinfo = ['> 900' if storage2DF.loc[storage2DF['BGC'] == iter1].Distance.tolist()[0] >= 900 else '< 900' for iter1 in superDict0[j].index.tolist()];
            
            # Removed dependency on storageDict
            bigSLICEinfo_v2 = ['> 900' if min(superDF2.loc[iter1]) >= 900 else '< 900' for iter1 in superDict0_j.index.tolist()]
            # find the distance to the closest GCF
            # closestGCFDist = [storage2DF.loc[storage2DF['BGC'] == iter1].Distance.tolist()[0] for iter1 in superDict0[j].index.tolist()];
            
            # Removed dependency on storageDict
            closestGCFDist_v2 = [min(superDF2.loc[iter1]) for iter1 in superDict0_j.index.tolist()]
            
            # Do some average distance calculations.
            for k in [1,2,3,4,5]:
                # find the distance to the X closest GCFs, then average them. Where X is set above as aveDist
                # aveDistStorage = [sorted(superDict0_j.loc[iter2])[:k] for iter2 in superDict0_j.index]; # list of lists
                # averageDist = [np.mean(entry1) for entry1 in aveDistStorage]; # average each list inside the list of lists
                # dfMDStemp['Average of ' + str(k) + ' Distances'] = averageDist;
                dfMDStemp[f'Average of {k} Distances'] = [np.mean(sorted(superDict0_j.loc[iter2])[:k]) for iter2 in superDict0_j.index]
    
            
            dfMDStemp['BIG-SLICE'] = bigSLICEinfo_v2; # not dependent on storageDict anymore
            dfMDStemp['Distance'] = closestGCFDist_v2; # not dependent on storageDict anymore
            # dfMDStemp['Average Distance'] = averageDist;
            # dfMDStemp['# of GCF Distances used'] = aveDist;
            # Add metadata
            # bgSCAPEinfo = [bgSCAPEclasses.loc[bgSCAPEclasses['name 1'] == i]['bgc product'].values for i in dfMDStemp.index];
            # dfMDStemp['Metadata'] = bgSCAPEinfo
            dfMDStemp['Metadata'] = dfMDStemp.index.map(lambda x: bgSCAPEclasses.loc[bgSCAPEclasses['name 1'] == x]['bgc product'].values[0])
            dfMDStemp.rename(columns={0:'x', 1:'y'}, inplace=True) # clean up column names for ease of reference later
            
            dictDFMDS[j][j+"_"+i] = dfMDStemp;
    
    specificMets = distMets # this is where the different strings for distance metrics go
    dictDFMDS_ALL['dictDFMDS_yesHybridCategory'] = dictDFMDS;
    # Save to Pickle (in current directory, which should be scriptOutputs)
    if not os.path.isfile('mds_data_run.pickle'):
        with open('mds_data_run.pickle', 'wb') as handle:
            pickle.dump(dictDFMDS_ALL, handle, protocol = pickle.HIGHEST_PROTOCOL)

# Timer
endTime = time.time()
lengthTime = endTime - startTime;
t_str = '{}h{}m{}s'.format(int(lengthTime/3600),int(lengthTime%3600/60),int(lengthTime%3600%60))
print('\nCompleted processing in {}'.format(t_str))

#%% Dash App
chdir(scriptOutputs)

colorSelectionScale = px.colors.named_colorscales()
classes1 = ['Terpene', 'PKS I', 'RiPP', 'PKS Other', 'PKS/NRPS', 'NRPS', 'Other', 'Saccharide', 'Hybrid (Not PKS/NRPS)']; # matches the colors indices.
substrings = ['Terpene', 'PKS I', 'RiPP', 'Other', 'PKS Other', 'PKS/NRPS', 'NRPS', 'Saccharide']
distMets = ['euclidean', 'cosine', 'cityblock', 'l1', 'l2', 'manhattan', 'braycurtis', 'canberra', 'chebyshev', 'correlation', 'hamming'] # this is where the different strings for distance metrics go


#%% Initialize App

# Sample empty plot
empty_plot = px.scatter()

app = dash.Dash(__name__)
app.title = "BGC Prioritization App"

# Define the layout of the app
app.layout = html.Div([
    # Left column with dropdowns and search bar
    html.Div([
        html.H1("BGC Prioritization Display", style={'textAlign': 'center'}),
        html.H3("Welcome to the Dash Tool to Assist with Prioritizing BGCs Based on Potential Novelty",
                style={'textAlign': 'center'}),
        html.H5("Developer: Imraan Alas (Bugni Lab)", style={'textAlign':'center'}),
        html.Div(id="intro",
                 children="Manipulate the various dropdown bars to customize which data is being presented, "
                          "the metadata used to color the plots, and the associated color scales.",
                 style={'textAlign': 'center'}),
        html.P("BGC Class: Choose from pre-defined BGC classes. You can't select Hybrid if Hybrids are not separated.", style={'textAlign': 'center'}),
        dcc.Dropdown(
            id='bigScapeClass-dropdown',
            options=[{'label': key, 'value': key} for key in classes1],
            value=substrings[0],  # Default value
            style={'width': '48%', 'margin': 'auto', 'text-align': 'center'},
            clearable=False,
            multi=False
        ),
        html.P("Distance Metric: Choose the distance metric used for metric dimensional scaling.", style={'textAlign': 'center'}),
        dcc.Dropdown(
            id='distanceMetric-dropdown',
            options=[{'label': metric, 'value': metric} for metric in distMets],
            value=distMets[0],  # Default value
            style={'width': '48%', 'margin': 'auto', 'text-align': 'center'},
            clearable=False,
            multi=False
        ),
        html.P("N Closest Distances: Choose the # of closest GCFs to a BGC that are considered for the average distance metadata.", style={'textAlign': 'center'}),
        dcc.Dropdown(
            id='averageDistance-dropdown',
            options=[{'label': metric, 'value': metric} for metric in [1, 2, 3, 4, 5]],
            value=3,  # Default value
            style={'width': '48%', 'margin': 'auto', 'text-align': 'center'},
            clearable=False,
            multi=False
        ),
        html.P("Color Palette: Choose a relevant color palette based on Plotly's inbuilt color scales.", style={'textAlign': 'center'}),
        dcc.Dropdown(
            id='colorPalette-dropdown',
            options=colorSelectionScale,
            value='rdpu',  # Default value
            style={'width': '48%', 'margin': 'auto', 'text-align': 'center'},
            clearable=False
        ),
        html.P("Hybrids Separated? Selecting 'Yes' will reduce the number of BGCs analyzed.", style={'textAlign':'center'}),
        dcc.Dropdown(
            id='hybridSelection-dropdown',
            options=[{'label': metric, 'value': metric} for metric in ['Yes', 'No']],
            value='No',  # Default value
            style={'width': '48%', 'margin': 'auto', 'text-align': 'center'},
            clearable=False
        ),
        html.P("Search BGC: Enter BGC name to filter the scatter plot.", style={'textAlign':'center'}),
        dcc.Input(
            id='search-input',
            type='text',
            value='',  # Default value
            placeholder='Enter BGC name...',
            style={'width': '99%', 'margin': 'auto', 'text-align': 'center'},
        ),
    ], style={'width': '48%', 'float': 'left'}),

    # Right column with the scatter plot
    html.Div([
        dash.html.Center(dcc.Graph(id='scatter-plot', figure=empty_plot)),
    ], style={'width': '48%', 'float': 'right'})

], style={"backgroundColor": "#f0f0f0", "display": "flex"})

# Define callback to update the scatter plot
@app.callback(
    Output('scatter-plot', 'figure'),
    [Input('bigScapeClass-dropdown', 'value'),
     Input('distanceMetric-dropdown', 'value'),
     Input('averageDistance-dropdown', 'value'),
     Input('colorPalette-dropdown', 'value'),
     Input('hybridSelection-dropdown', 'value'),
     Input('search-input', 'value')]
)
def update_scatter_plot(selected_big_scape_class, selected_distance_metric, selected_average_distance,
                         selected_color_scale, selected_hybrids, search_value):
    try:
        selected_key = f"{selected_big_scape_class}_orig_{selected_distance_metric}"
        if selected_hybrids == 'Yes':
            selectedHybridVal = 'yes'
        else:
            selectedHybridVal = 'no'
    except:
        return empty_plot

    try:
        if selected_key in dictDFMDS_ALL['dictDFMDS_' + selectedHybridVal + 'HybridCategory'][selected_big_scape_class + '_orig']:
            bgcTitleName = selected_key.split('_')[0]
            metricTitleName = selected_key.split('_')[-1].capitalize()
            df = dictDFMDS_ALL['dictDFMDS_' + selectedHybridVal + 'HybridCategory'][selected_big_scape_class + '_orig'][selected_key]

            selected_aveDist = 'Average of ' + str(selected_average_distance) + ' Distances'
            fig = px.scatter(df, x='x', y='y', color_continuous_scale=selected_color_scale, color=selected_aveDist,
                             size=selected_aveDist, hover_data=['Metadata', 'Distance'],
                             hover_name=['WMM' + i.capitalize() for i in df.index.tolist()])
            fig.update_layout(height=900, width=1000)
            fig.update_layout(title_text=f'{bgcTitleName} BGCs Analyzed with {metricTitleName} Distance Metric '
                                          f'<br><sup>Size & Color Described by Average Distance of BGC to Closest {selected_average_distance} GCFs</sup>',
                              title_x=0.5)
            fig.update_layout(xaxis_title='Dimension 1', yaxis_title='Dimension 2')

            if search_value:
                df['IsSearchMatch'] = df.index.str.contains(search_value, case=False)
                fig.update_traces(
                    marker=dict(line=dict(color='black', width=[3 if is_search_match and search_value else 0.5 for is_search_match in df['IsSearchMatch']])))
            else:
                fig.update_traces(marker=dict(line=dict(color='white', width=1)))
            return fig
        else:
            return empty_plot
    except:
        return empty_plot

# Run the app
if __name__ == '__main__':
    print('Open localhost:8002 or 127.0.0.1:8002 to see the dash visualization')
    app.run_server(port=8002)  # Open 127.0.0.1:8002 or localhost:8002 on your local web browser.
