

from multiprocessing import pool
from os import listdir
from os.path import isfile, join
from collections import defaultdict
import csv
import os
import io
from pprint import pprint
import math
from sqlite3 import Row
import numpy as np # python3 -m pip install -U pip
import pandas as pd
from io import StringIO
import pandas as pd



#######################################################################
"""STANDARDS"""
#######################################################################

# create dataframe with all standards

standard_df = []
standard_well = ["A","B","C"]
for val in range(1,7):
    temp_df = []
    for i in standard_well:
        standard1 = "{}{}".format(i,val)
        temp_df.append(standard1)
    standard_df.append(temp_df)
        
#print(standard_df)
dict_standard = defaultdict()

path = "/Users/aotgonbayar/Desktop/Python/"
os.chdir(path)
        
key = 0
cq_current = []
cq_list = []

for list_df in standard_df:

    key+=1

    with open('qPCR_WellResult_Run2.csv', encoding='unicode_escape') as file:
        data = csv.DictReader(file)
    
        for row in data: # each row in data
            if "Std{}".format(key) in row["Sample"]:
                
            
                dict_standard[key] = []
                cq_current = row["Cq"]
                cq_list.append(cq_current)
        dict_standard[key].append(cq_list)
        cq_list = []
#print(dict_standard) #has all the standards 1,2,3,4,5,6

# Replicate data points should differ by <=-0.2 cycles. Pull out outliers otherwise or repeat experiment:
#  
cq_delta_standard= defaultdict()
for standard_key,standard_list in dict_standard.items():
    #print(standard_key,standard_list)
    #print(standard_list[0])
    cq1 = float(standard_list[0][0]) # first value in each standard list
    cq2 = float(standard_list[0][1]) # second value in each standard list
    cq3 = float(standard_list[0][2]) # third value in each standard list
    cq_12 = abs(cq1-cq2)
    cq_13 = abs(cq1-cq3)
    cq_23 = abs(cq2-cq3)
    cq_compiled = [cq_12, cq_13, cq_23]
    cq_delta_standard[standard_key] = []
    cq_delta_standard[standard_key].append(cq_compiled)
    cq_compiled = []
#print(cq_delta_standard) # all the delta cq for all standard key

dict_standard_no_out = defaultdict()
for standard_key, st_list in cq_delta_standard.items():
    #print(standard_key,st_list)
    key_list = cq_delta_standard[standard_key][0] # gets the individual list per key
    outlier_idx = []
    for idx, val in enumerate(key_list): # access individual value per list
        if val>0.2:
            outlier_idx += [idx]
            key_list[idx] = 0
            """outlier = val
            find_true_outlier = []
            find_true_outlier.append(outlier) # 
    if outlier in key_list: # if the outlier is in the list, we want to remove the original standard
        idx = key_list.index(outlier) # find the index of this outlier from the delta cq list
        key_list[idx] = 0 # replace value with 0 to keep the same index"""
            idx_new_max = key_list.index(max(key_list))
            outlier_idx += [idx_new_max]
            key_list[idx_new_max] = 0 # replace value with 0 to keep same index
    #print(outlier_idx)
    if not outlier_idx:
        dict_standard_no_out[standard_key] = dict_standard[standard_key][0]
    else:
        if set(outlier_idx) == {0,1}:
            outlier_idx = 0
        elif set(outlier_idx) == {0,2}:
            outlier_idx = 1
        elif set(outlier_idx) == {1,2}:
            outlier_idx = 2
        else:
            print(f"ERROR Multiple outliers for standard {standard_key}")
        #print(outlier_idx)
        new_list = dict_standard[standard_key][0].copy()
        #print(new_list)
        del new_list[outlier_idx]
        print(new_list)
        dict_standard_no_out[standard_key] = new_list
#print(dict_standard_no_out)

    
        # Recall cq_compiled = [cq_12, cq_13, cq_23]
        # Recall key_list = [0,1,2]
    # if there is 0 on the following indices of key_list:
        # key_list = [0,0,x] --> corresponds w/ cq_12, cq_13; thus cq1 is the outlier
        # key_list = [0,x,0] --> corresponds w/ cq_12, cq_23; thus cq2 is the outlier
        # key_list = [x,0,0] --> corresponds w/ cq_13, cq_23; thus cq3 is the outlier
        # Note: IF there is more than one outlier from the triplicate; print out "ERROR message"
    # now, we want to recall this:
    # cq1 = float(standard_list[0][0]) # first value in each standard list
    # cq2 = float(standard_list[0][1]) # second value in each standard list
    # cq3 = float(standard_list[0][2]) # third value in each standard list
        # this will help locate the outlier


#######################################################################
"""SAMPLES"""
#######################################################################

standard_df = []
standard_well = ["A","B","C"]
for val in range(1,7):
    temp_df = []
    for i in standard_well:
        standard1 = "{}{}".format(i,val)
        temp_df.append(standard1)
    standard_df.append(temp_df)
        
#print(standard_df)
dict_standard = defaultdict()

path = "/Users/aotgonbayar/Desktop/Python/"
os.chdir(path)
        
key = 0
cq_current = []
cq_list = []

for list_df in standard_df:

    key+=1

    with open('qPCR_WellResult_Run2.csv', encoding='unicode_escape') as file:
        data = csv.DictReader(file)
    
        for row in data: # each row in data
            if "Std{}".format(key) in row["Sample"]:
                
            
                dict_standard[key] = []
                cq_current = row["Cq"]
                cq_list.append(cq_current)
        dict_standard[key].append(cq_list)
        cq_list = []
#print(dict_standard) #has all the standards 1,2,3,4,5,6

# Replicate data points should differ by <=-0.2 cycles. Pull out outliers otherwise or repeat experiment:
#  
cq_delta_standard= defaultdict()
for standard_key,standard_list in dict_standard.items():
    #print(standard_key,standard_list)
    #print(standard_list[0])
    cq1 = float(standard_list[0][0]) # first value in each standard list
    cq2 = float(standard_list[0][1]) # second value in each standard list
    cq3 = float(standard_list[0][2]) # third value in each standard list
    cq_12 = abs(cq1-cq2)
    cq_13 = abs(cq1-cq3)
    cq_23 = abs(cq2-cq3)
    cq_compiled = [cq_12, cq_13, cq_23]
    cq_delta_standard[standard_key] = []
    cq_delta_standard[standard_key].append(cq_compiled)
    cq_compiled = []
#print(cq_delta_standard) # all the delta cq for all standard key

dict_standard_no_out = defaultdict()
for standard_key, st_list in cq_delta_standard.items():
    #print(standard_key,st_list)
    key_list = cq_delta_standard[standard_key][0] # gets the individual list per key
    outlier_idx = []
    for idx, val in enumerate(key_list): # access individual value per list
        if val>0.2:
            outlier_idx += [idx]
            key_list[idx] = 0
            """outlier = val
            find_true_outlier = []
            find_true_outlier.append(outlier) # 
    if outlier in key_list: # if the outlier is in the list, we want to remove the original standard
        idx = key_list.index(outlier) # find the index of this outlier from the delta cq list
        key_list[idx] = 0 # replace value with 0 to keep the same index"""
            idx_new_max = key_list.index(max(key_list))
            outlier_idx += [idx_new_max]
            key_list[idx_new_max] = 0 # replace value with 0 to keep same index
    #print(outlier_idx)
    if not outlier_idx:
        dict_standard_no_out[standard_key] = dict_standard[standard_key][0]
    else:
        if set(outlier_idx) == {0,1}:
            outlier_idx = 0
        elif set(outlier_idx) == {0,2}:
            outlier_idx = 1
        elif set(outlier_idx) == {1,2}:
            outlier_idx = 2
        else:
            print(f"ERROR Multiple outliers for standard {standard_key}")
        #print(outlier_idx)
        new_list = dict_standard[standard_key][0].copy()
        #print(new_list)
        del new_list[outlier_idx]
        print(new_list)
        dict_standard_no_out[standard_key] = new_list
print(dict_standard_no_out)

