import sys
import os
from collections import Counter

# Parse a .fasta file
def parse_fasta(fname):
    identifiers = []
    with open(fname, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                identifiers.append(line[1:].strip())
    return identifiers
    
# For the given path, get the List of all files in the directory tree 
def getListOfFiles(dirName):
    # create a list of file and sub directories 
    # names in the given directory 
    listOfFile = os.listdir(dirName)
    allFiles = list()
    # Iterate over all the entries
    for entry in listOfFile:
        # Create full path
        fullPath = os.path.join(dirName, entry)
        # If entry is a directory then get the list of files in this directory 
        if os.path.isdir(fullPath):
            allFiles = allFiles + getListOfFiles(fullPath)
        else:
            allFiles.append(fullPath)
                
    return allFiles

def give_dataset_info(metadata_df, voc):

    pango_lin = 'pangolin_lineage'

    if pango_lin not in metadata_df:
        pango_lin = "pango_lineage"

    # unique lineages
    pango_lin_categories = set(metadata_df[pango_lin])
    # amount per unique lineage
    counts_per_pango_lin = Counter(metadata_df[pango_lin])
    # total amount of sequences
    sum_counts = sum(counts_per_pango_lin.values())
    # Categories of lineages found
    unique_lineages_found_amount =  len(pango_lin_categories)
    # avg variation per lineage
    avg_variation_per_lineage = round(sum_counts/ len(pango_lin_categories), 2)
    # counts for voc
    counts_for_voc = Counter(metadata_df[pango_lin]).get(voc)
    # find min max dates
    max_date = max(metadata_df['date'])
    min_date = min(metadata_df['date'])

    return unique_lineages_found_amount, avg_variation_per_lineage, counts_for_voc,max_date, min_date, pango_lin_categories 

def info_to_dict(info_df):
    info_df['amount'] = [float(x) for x in info_df['amount']]
    keys = set(info_df['clade'])
    dict_clade_counts  = dict.fromkeys(keys, 0)

    total = 0
    for row in info_df.itertuples():
        key = row[3]
        if key == "None":
            continue
        count = row[2]
        existing_cnt = dict_clade_counts.get(key)
        new_cnt = existing_cnt + count
        dict_clade_counts[key] = new_cnt
        total = row[4]

    print(dict_clade_counts)
    print(sum(dict_clade_counts.values()),total )

    return dict_clade_counts