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
import statistics


# samples !!
standard_df = []
standard_well = ["D","E","F"] # triplicates
for val in range(1,10): # 9 samples , val = 1-9 inclusive
    temp_df = []
    for i in standard_well: #i =  [D, E F] x9 samples
        standard1 = "{}{}".format(i,val)
        temp_df.append(standard1)
    standard_df.append(temp_df)
        
#print(standard_df) # all sample wells by name
dict_standard = defaultdict()
dict_qty = defaultdict()

path = "/Users/aotgonbayar/Desktop/Python/"
os.chdir(path)
        
key = 0
cq_current = []
cq_list = []
qty_list = []
qty_current = []

for list_df in standard_df: #[D1,E1,F1], etc
    
    
    key+=1
    

    with open('qPCR_WellResult_Run2.csv', encoding='unicode_escape') as file:
        data = csv.DictReader(file)
    
        for row in data: # each row in data
            if list_df[0] in row["Well"]: # first value

                dict_standard[key] = []
                dict_qty[key] = []
                cq_current = row["Cq"]
                qty_current = row["Quantity"]

                cq_list.append(cq_current)
                qty_list.append(qty_current)
        #dict_standard[key].append(cq_list)
        #cq_list = []
            elif list_df[1] in row["Well"]: # second value


                cq_current = row["Cq"]
                qty_current = row["Quantity"]

                cq_list.append(cq_current)
                qty_list.append(qty_current)

            elif list_df[2] in row["Well"]: # second value

                cq_current = row["Cq"]
                qty_current = row["Quantity"]

                cq_list.append(cq_current)
                qty_list.append(qty_current)

        dict_standard[key].append(cq_list)
        dict_qty[key].append(qty_list)

        
        cq_list = []
        qty_list = []


#print(dict_standard) #has all the samples 1,2,3,4,5,6,7,8,9 relative to plate
#print(dict_qty)

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
dict_qty_no_out = defaultdict()
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
        dict_qty_no_out[standard_key] = dict_qty[standard_key][0]
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
        new_list2 = dict_qty[standard_key][0].copy()
        #print(new_list)
        del new_list[outlier_idx]
        del new_list2[outlier_idx]
        #print(new_list)
        dict_standard_no_out[standard_key] = new_list
        dict_qty_no_out[standard_key] = new_list2
#print(dict_standard_no_out) # has all original cq values without the outliers
#print(dict_qty_no_out)
average_pM = []

for key, qty_list in dict_qty_no_out.items():
    qty_list_float = [float(n) for n in qty_list]
    
    average_qty_current_key = statist=statistics.mean((qty_list_float))
    print(average_qty_current_key)

#df = pd.DataFrame(dict_standard)
#print((df))

