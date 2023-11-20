# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 21:22:06 2023

@author: imraan
"""

# NOTE: This file contains scripts on how to scrape the 30S and 50S ribosomal protein nucleotides and amino acids, in compatible form for MaffT.

#%%

# original packages:
import os
from bs4 import BeautifulSoup # beautiful soup version 4.10.0
from copy import deepcopy 
from os import chdir, getcwd
import re

import pandas as pd # pandas version 1.3.5
import numpy as np # numpy version 1.21.6
import matplotlib.pyplot as plt
import json

#%% 

chdir('C:/Users/imraa/Documents/UWM/BugniLab/Genome Assembly/bakta/baktaOutputs/20231102-murphyLab')    

#%% Nucleotide Version:

# Find the files.
def find_json_files(directory = '.'):
    json_files = [];
    # Walk through the directory and its subdirectories
    for root, dirs, files in os.walk(directory):
        # Check if there are any JSON files in the current directory
        json_files.extend([os.path.join(root, file) for file in files if file.endswith('.json')])
    return json_files

json_files_list = find_json_files()

# Open the text file that everything is going to go into:
with open("nucleotide_ribosomalProteins.fasta","w") as megaFile:
    # Choose the json files!
    for json_file in json_files_list:
        genomeName = json_file.split('.json')[0].split('\\')[-1] # hard-coded.
        # Open the file and load the JSON data
        with open(json_file, 'r') as file:
            data = json.load(file)
        # Create the list.
        filtered_data = [];
        # Scrape the relevant information. (50S or 30S)
        for key, item in data.items():
            if isinstance(item, list):
                for index, feat in enumerate(item):
                    if isinstance(feat, dict) and feat.get("product", "").startswith(("30S", "50S")):
                        index_in_data = data['features'].index(feat)
                        filtered_data.append({
                            "key": key,
                            "index_in_data": index_in_data,
                            "feature_data": feat,
                            # Include other relevant keys/values as needed
                        })
        # remove the 50S or 30S that isn't the ribosomal protein (example: methyltransferases)
        # Simple filter: Gotta have a number at the end of it.
        filtered_data_pure = [item for item in filtered_data if item['feature_data']["product"].endswith(tuple("0123456789"))];
        # Write everything to the file.
        for item in filtered_data_pure:
            rel_info = item['feature_data'];
            symbol = ">";
            product = rel_info.get("product","").replace(" ", "_");
            contig = rel_info.get("contig","")
            start = rel_info.get("start","")
            stop = rel_info.get("stop","")
            strand = rel_info.get("strand","")
            nt = rel_info.get("nt","")
            # print(f"{symbol}{genomeName}_{product}::{contig}:{start}-{stop}({strand})\n")
            megaFile.write(f"{symbol}{genomeName}_{product}::{contig}:{start}-{stop}({strand})\n")
            megaFile.write(f"{nt}\n")


#%% Aminoacid Version:

# Find the files.
def find_json_files(directory = '.'):
    json_files = [];
    # Walk through the directory and its subdirectories
    for root, dirs, files in os.walk(directory):
        # Check if there are any JSON files in the current directory
        json_files.extend([os.path.join(root, file) for file in files if file.endswith('.json')])
    return json_files

json_files_list = find_json_files()

# Open the text file that everything is going to go into:
with open("aminoAcid_ribosomalProteins.fasta","w") as megaFile: # CHANGE FILE NAME
    # Choose the json files!
    for json_file in json_files_list:
        genomeName = json_file.split('.json')[0].split('\\')[-1] # hard-coded.
        # Open the file and load the JSON data
        with open(json_file, 'r') as file:
            data = json.load(file)
        # Create the list.
        filtered_data = [];
        # Scrape the relevant information. (50S or 30S)
        for key, item in data.items():
            if isinstance(item, list):
                for index, feat in enumerate(item):
                    if isinstance(feat, dict) and feat.get("product", "").startswith(("30S", "50S")):
                        index_in_data = data['features'].index(feat)
                        filtered_data.append({
                            "key": key,
                            "index_in_data": index_in_data,
                            "feature_data": feat,
                            # Include other relevant keys/values as needed
                        })
        # remove the 50S or 30S that isn't the ribosomal protein (example: methyltransferases)
        # Simple filter: Gotta have a number at the end of it.
        filtered_data_pure = [item for item in filtered_data if item['feature_data']["product"].endswith(tuple("0123456789"))];
        # Write everything to the file.
        for item in filtered_data_pure:
            rel_info = item['feature_data'];
            symbol = ">";
            product = rel_info.get("product","").replace(" ", "_");
            contig = rel_info.get("contig","")
            start = rel_info.get("start","")
            stop = rel_info.get("stop","")
            strand = rel_info.get("strand","")
            aa = rel_info.get("aa","")
            # print(f"{symbol}{genomeName}_{product}::{contig}:{start}-{stop}({strand})\n")
            megaFile.write(f"{symbol}{genomeName}_{product}::{contig}:{start}-{stop}({strand})\n")
            megaFile.write(f"{aa}\n")