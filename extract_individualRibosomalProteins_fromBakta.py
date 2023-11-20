# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 21:22:06 2023

@author: imraan
"""

# NOTE: This file contains scripts on how to scrape the 30S and 50S ribosomal protein nucleotides and amino acids, in compatible form for MaffT.
# The previous version of this lumped everything into one file. This one generates an individual file for each protein.

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
def find_json_files(directory='.'):
    json_files = []
    # Walk through the directory and its subdirectories
    for root, dirs, files in os.walk(directory):
        # Check if there are any JSON files in the current directory
        json_files.extend([os.path.join(root, file) for file in files if file.endswith('.json')])
    return json_files

json_files_list = find_json_files()

# Dictionary to store sequences for each protein
protein_sequences = {}

# Iterate through the JSON files
for json_file in json_files_list:
    genome_name = json_file.split('.json')[0].split('\\')[-1]  # hard-coded.
    # Open the file and load the JSON data
    with open(json_file, 'r') as file:
        data = json.load(file)
    # Create the list.
    filtered_data = []
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
    # Remove the 50S or 30S that isn't the ribosomal protein (example: methyltransferases)
    # Simple filter: Gotta have a number at the end of it.
    filtered_data_pure = [item for item in filtered_data if item['feature_data']["product"].endswith(tuple("0123456789"))]

    # Process the filtered data
    for item in filtered_data_pure:
        rel_info = item['feature_data']
        product = rel_info.get("product", "").replace(" ", "_")

        # Append a small identifying element to the protein name if already present
        if product not in protein_sequences:
            protein_sequences[product] = []

        nt = rel_info.get("nt", "")
        # Append the sequence to the list for that protein
        protein_sequences[product].append(
            f">{genome_name}_{product}::{rel_info.get('contig', '')}:{rel_info.get('start', '')}-{rel_info.get('stop', '')}({rel_info.get('strand', '')})")
        protein_sequences[product].append(nt)

# Iterate through the protein_sequences dictionary
for protein, sequences in protein_sequences.items():
    # Replace problematic characters in the protein name for the file
    safe_protein_name = protein.replace("/", "_")  # L7/L12 is not a thing Windows likes as a file name.

    # Create the output file name
    output_file_name = f"20231102_nucleotide_{safe_protein_name}.fasta"

    # Write sequences to the output file
    with open(output_file_name, "w") as output_file:
        output_file.write("\n".join(sequences))

#%% Aminoacid Version:

# Find the files.
def find_json_files(directory='.'):
    json_files = []
    # Walk through the directory and its subdirectories
    for root, dirs, files in os.walk(directory):
        # Check if there are any JSON files in the current directory
        json_files.extend([os.path.join(root, file) for file in files if file.endswith('.json')])
    return json_files

json_files_list = find_json_files()

# Dictionary to store sequences for each protein
protein_sequences = {}

# Iterate through the JSON files
for json_file in json_files_list:
    genome_name = json_file.split('.json')[0].split('\\')[-1]  # hard-coded.
    # Open the file and load the JSON data
    with open(json_file, 'r') as file:
        data = json.load(file)
    # Create the list.
    filtered_data = []
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
    # Remove the 50S or 30S that isn't the ribosomal protein (example: methyltransferases)
    # Simple filter: Gotta have a number at the end of it.
    filtered_data_pure = [item for item in filtered_data if item['feature_data']["product"].endswith(tuple("0123456789"))]

    # Process the filtered data
    for item in filtered_data_pure:
        rel_info = item['feature_data']
        product = rel_info.get("product", "").replace(" ", "_")

        # Append a small identifying element to the protein name if already present
        if product not in protein_sequences:
            protein_sequences[product] = []

        aa = rel_info.get("aa", "")
        # Append the sequence to the list for that protein
        protein_sequences[product].append(
            f">{genome_name}_{product}::{rel_info.get('contig', '')}:{rel_info.get('start', '')}-{rel_info.get('stop', '')}({rel_info.get('strand', '')})")
        protein_sequences[product].append(aa)

# Iterate through the protein_sequences dictionary
for protein, sequences in protein_sequences.items():
    # Replace problematic characters in the protein name for the file
    safe_protein_name = protein.replace("/", "_")  # L7/L12 is not a thing Windows likes as a file name.

    # Create the output file name
    output_file_name = f"20231102_aminoacid_{safe_protein_name}.fasta"

    # Write sequences to the output file
    with open(output_file_name, "w") as output_file:
        output_file.write("\n".join(sequences))

