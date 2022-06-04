# Goal: Determine for potential overloading of cells due to too many PCR cycles

# initialize packages 
from os import listdir
from os.path import isfile, join
from collections import defaultdict
import csv
import os
import io
from pprint import pprint
import math
import numpy as np # python3 -m pip install -U pip
import pandas as pd

# 1. Organize Data by Pool 1, 2, 3, 4
  # Pool:[Size (bp), Concentration (ng/ul), Markers]

#def hashmap_per_batch(start_file_name):
path = "/Users/aotgonbayar/Desktop/Python/"
a  = os.listdir(path)  # lists all files in this folder
from io import StringIO

# dir_path = os.path.dirname(os.path.realpath(__file__)) To get the full path to the directory a Python file is contained in, write this in that file:



hash_map_pool_1 = defaultdict(list)
hash_map_pool_2 = defaultdict(list)
hash_map_pool_3 = defaultdict(list)
hash_map_pool_4 = defaultdict(list)


for csv_doc in a:# list all files in selected folder
    current_file = os.listdir(path)
    
    
    #print(type(current_file[0]))
    if csv_doc.__contains__('Batch38'): # if csv file contains Peak Table
        #print(csv_doc)
        #print(current_file)
        #print(os.getcwd())
        #print(csv_doc)        #print(trial_number)
        os.chdir(path) # change working directory
        #print(os.getcwd)
        #with open(current_file[0],"rt", encoding = "utf8") as file: 
        #with open(csv_doc,"rt", encoding = "utf8") as file:
          #data = csv.DictReader(file)
          #data = pd.read_csv(current_file[0], encoding = 'unicode_escape')
        
        
        data = pd.read_csv(csv_doc, encoding = 'unicode_escape') 
        #print(os.getcwd())

        # look through data and find filter out rows with ladder
        for count, row in enumerate(data.Well): # All rows under "Well" data frame
          #print(row)
        
        #HashMap[""][trial_number] = []
          #Separate all pools
          if row=="EL1": # ignore electronic ladder data
              continue

          elif "1" in row: # A1,B1,C1
            pool_data = data.iloc[count] # pull out all the rows
            size = pool_data[3] # units in bp
            concentration = pool_data[4] # units in ng/ul

            markers = pool_data.Observations # peak information

            
            # add to dictionary
            hash_map_pool_1["Size"].append(size)
            hash_map_pool_1["Concentration"].append(concentration)
            hash_map_pool_1["Markers"].append(markers)
            # POOL 1 DATA: hash_map_pool_1

          elif "2" in row: # A2, B2, C2
            pool_data = data.iloc[count] # pull out all the rows
            size = pool_data[3] # units in bp
            concentration = pool_data[4] # units in ng/ul
            markers = pool_data.Observations # peak information
            # add to dictionary
            hash_map_pool_2["Size"].append(size)
            hash_map_pool_2["Concentration"].append(concentration)
            hash_map_pool_2["Markers"].append(markers)
            # POOL 2 DATA: hash_map_pool_2

          elif "3" in row: # A3, B3, C3
            pool_data = data.iloc[count] # pull out all the rows
            size = pool_data[3] # units in bp
            concentration = pool_data[4] # units in ng/ul
            markers = pool_data.Observations # peak information
            # add to dictionary
            hash_map_pool_3["Size"].append(size)
            hash_map_pool_3["Concentration"].append(concentration)
            hash_map_pool_3["Markers"].append(markers)
            # POOL 3 DATA: hash_map_pool_3


          elif "4" in row: # A4, B4, C4
            pool_data = data.iloc[count] # pull out all the rows
            size = pool_data[3] # units in bp
            concentration = pool_data[4] # units in ng/ul
            markers = pool_data.Observations # peak information
            # add to dictionary
            hash_map_pool_4["Size"].append(size)
            hash_map_pool_4["Concentration"].append(concentration)
            hash_map_pool_4["Markers"].append(markers)
  
              # POOL 4 DATA: hash_map_pool_4
              
#print(hash_map_pool_1)
#print(hash_map_pool_2)
#print(hash_map_pool_3)
#print(hash_map_pool_4)
                                          #######################################
#2. Organize hashmaps by peak
# Find value between Upper marker and Lower Marker 



# FINAL VALUES
pool_1 = defaultdict(list)
pool_2 = defaultdict(list)
pool_3 = defaultdict(list)
pool_4 = defaultdict(list)
#print(hash_map_pool_1["Markers"])
#############################################POOL 1####################################################################

counter = len(hash_map_pool_1["Markers"]) # total number of values in list, of dictionary
#print(counter)
idx_remove = [] # create empty list

for tracker, val in enumerate(hash_map_pool_1["Markers"]): 
    val = str(val)
    if " Marker" in val:
    # find the index of this occurence
        idx_remove.append(tracker)
        #print(tracker) # idx_remove has all values "lower to upper marker" #
        #  Now get individually calculated peaks. len(idx_remove) should always be even!

#print(idx_remove)
peak1 = [] # row with lower marker and upper marker to find values
find_peak1 = []  # ALL PEAK VALUES !! per pool
peak_val1 = []
idx_dict = []
for r in idx_remove:
  #print(r)
  if "Lower Marker" in hash_map_pool_1["Markers"][r]: # instances with Lower Marker
    peak1.append(r) # append lower marker index
    # now find upper marker index 
  if "Upper Marker" in hash_map_pool_1["Markers"][r]: # instances with Upper Marker
    peak1.append(r) # append upper marker index
    #if len(peak) ==2:
    #print(peak)
    find_peak1 = list(range(peak1[0],peak1[1]))
    find_peak1.pop(0)
    #print(find_peak1) # All peak idx for this pool
    idx_dict.append(find_peak1)
    peak1 = [] # start new
    find_peak1 = []

#print(idx_dict)  # all index list

# pull out values using index from idx_dict from hash_map to find peaks
pks = []

for key_index in idx_dict: # key_index: [1,2],[5,6]

  for tracker_key, key in enumerate(hash_map_pool_1["Concentration"]): #key are all values
    #print(tracker_key,key) # tracker is the index of hashmap
    if tracker_key in key_index: # if the index is in the concentration, pull it out
      pks.append(key)
  pool_1["Concentration_Peak1"].append(pks)
  pks = []

  for tracker_key_1, key_1 in enumerate(hash_map_pool_1["Size"]):
    if tracker_key_1 in key_index: # if the index is in the concentration, pull it out
      pks.append(key_1)

  pool_1["Size"].append(pks)
  pks = []
#print(pool_1) #has all peak values from pool 1 hashmap
###############################################POOL 2############################################################################

counter = len(hash_map_pool_2["Markers"]) # total number of values in list, of dictionary
#print(counter)
idx_remove = [] # create empty list

for tracker, val in enumerate(hash_map_pool_2["Markers"]): 
    val = str(val)
    if " Marker" in val:
    # find the index of this occurence
        idx_remove.append(tracker)
        #print(tracker) # idx_remove has all values "lower to upper marker" #
        #  Now get individually calculated peaks. len(idx_remove) should always be even!

#print(idx_remove)
peak1 = [] # row with lower marker and upper marker to find values
find_peak1 = []  # ALL PEAK VALUES !! per pool
peak_val1 = []
idx_dict = []
for r in idx_remove:
  #print(r)
  if "Lower Marker" in hash_map_pool_2["Markers"][r]: # instances with Lower Marker
    peak1.append(r) # append lower marker index
    # now find upper marker index 
  if "Upper Marker" in hash_map_pool_2["Markers"][r]: # instances with Upper Marker
    peak1.append(r) # append upper marker index
    #if len(peak) ==2:
    #print(peak)
    find_peak1 = list(range(peak1[0],peak1[1]))
    find_peak1.pop(0)
    #print(find_peak1) # All peak idx for this pool
    idx_dict.append(find_peak1)
    peak1 = [] # start new
    find_peak1 = []

#print(idx_dict)  # all index list

# pull out values using index from idx_dict from hash_map to find peaks
pks = []

for key_index in idx_dict: # key_index: [1,2],[5,6]

  for tracker_key, key in enumerate(hash_map_pool_2["Concentration"]): #key are all values
    #print(tracker_key,key) # tracker is the index of hashmap
    if tracker_key in key_index: # if the index is in the concentration, pull it out
      pks.append(key)
  pool_2["Concentration_Peak1"].append(pks)
  pks = []

  for tracker_key_1, key_1 in enumerate(hash_map_pool_2["Size"]):
    if tracker_key_1 in key_index: # if the index is in the concentration, pull it out
      pks.append(key_1)

  pool_2["Size"].append(pks)
  pks = []
#print(pool_2) #has all peak values from pool 1 hashmap

###############################################POOL 3############################################################################




counter = len(hash_map_pool_3["Markers"]) # total number of values in list, of dictionary
#print(counter)
idx_remove = [] # create empty list

for tracker, val in enumerate(hash_map_pool_3["Markers"]): 
    val = str(val)
    if " Marker" in val:
    # find the index of this occurence
        idx_remove.append(tracker)
        #print(tracker) # idx_remove has all values "lower to upper marker" #
        #  Now get individually calculated peaks. len(idx_remove) should always be even!

#print(idx_remove)
peak1 = [] # row with lower marker and upper marker to find values
find_peak1 = []  # ALL PEAK VALUES !! per pool
peak_val1 = []
idx_dict = []
for r in idx_remove:
  #print(r)
  if "Lower Marker" in hash_map_pool_3["Markers"][r]: # instances with Lower Marker
    peak1.append(r) # append lower marker index
    # now find upper marker index 
  if "Upper Marker" in hash_map_pool_3["Markers"][r]: # instances with Upper Marker
    peak1.append(r) # append upper marker index
    #if len(peak) ==2:
    #print(peak)
    find_peak1 = list(range(peak1[0],peak1[1]))
    find_peak1.pop(0)
    #print(find_peak1) # All peak idx for this pool
    idx_dict.append(find_peak1)
    peak1 = [] # start new
    find_peak1 = []

#print(idx_dict)  # all index list

# pull out values using index from idx_dict from hash_map to find peaks
pks = []

for key_index in idx_dict: # key_index: [1,2],[5,6]

  for tracker_key, key in enumerate(hash_map_pool_3["Concentration"]): #key are all values
    #print(tracker_key,key) # tracker is the index of hashmap
    if tracker_key in key_index: # if the index is in the concentration, pull it out
      pks.append(key)
  pool_3["Concentration_Peak1"].append(pks)
  pks = []

  for tracker_key_1, key_1 in enumerate(hash_map_pool_3["Size"]):
    if tracker_key_1 in key_index: # if the index is in the concentration, pull it out
      pks.append(key_1)

  pool_3["Size"].append(pks)
  pks = []
#print(pool_3) #has all peak values from pool 3 hashmap

###############################################POOL 4############################################################################

counter = len(hash_map_pool_4["Markers"]) # total number of values in list, of dictionary
#print(counter)
idx_remove = [] # create empty list

for tracker, val in enumerate(hash_map_pool_4["Markers"]): 
    val = str(val)
    if " Marker" in val:
    # find the index of this occurence
        idx_remove.append(tracker)
        #print(tracker) # idx_remove has all values "lower to upper marker" #
        #  Now get individually calculated peaks. len(idx_remove) should always be even!

#print(idx_remove)
peak1 = [] # row with lower marker and upper marker to find values
find_peak1 = []  # ALL PEAK VALUES !! per pool
peak_val1 = []
idx_dict = []
for r in idx_remove:
  #print(r)
  if "Lower Marker" in hash_map_pool_4["Markers"][r]: # instances with Lower Marker
    peak1.append(r) # append lower marker index
    # now find upper marker index 
  if "Upper Marker" in hash_map_pool_4["Markers"][r]: # instances with Upper Marker
    peak1.append(r) # append upper marker index
    #if len(peak) ==2:
    #print(peak)
    find_peak1 = list(range(peak1[0],peak1[1]))
    find_peak1.pop(0)
    #print(find_peak1) # All peak idx for this pool
    idx_dict.append(find_peak1)
    peak1 = [] # start new
    find_peak1 = []

#print(idx_dict)  # all index list

# pull out values using index from idx_dict from hash_map to find peaks
pks = []

for key_index in idx_dict: # key_index: [1,2],[5,6]

  for tracker_key, key in enumerate(hash_map_pool_4["Concentration"]): #key are all values
    #print(tracker_key,key) # tracker is the index of hashmap
    if tracker_key in key_index: # if the index is in the concentration, pull it out
      pks.append(key)
  pool_4["Concentration_Peak1"].append(pks)
  pks = []

  for tracker_key_1, key_1 in enumerate(hash_map_pool_4["Size"]):
    if tracker_key_1 in key_index: # if the index is in the concentration, pull it out
      pks.append(key_1)

  pool_4["Size"].append(pks)
  pks = []
#print(pool_4) #has all peak values from pool 4

#### summary
print(pool_1)
print(pool_2)
print(pool_3)
print(pool_4)


#STEP 3: Write data from pool into CSV file for analysis

# open csv file for writing
# f = open('/Users/aotgonbayar/Desktop/Python/code.csv','w') # open the file in the write mode

# # create the csv writer
# writer_csv = csv.writer(f)

# # write a row to the csv file
# writer_csv.writerow(row)

# # close the file
# f.close()


  

#for indx, xyz in enumerate(hash_map_pool_1["Concentration"]):
  # if idx_dict index is not in this dictionary, remove it
  #print(indx, xyz) # indx is the index [0,1,2,3,4,5,6,7,8,9,10,11], xyz is the value of the index [7.5,1.3.....]
  #if indx not in idx_dict["Index"][0] or indx not in idx_dict["Index"][1]
  

  #if idx_dict["Index"][indx] not in tuple(idx_dict):
    #hash_map_pool_1["Concentration"].pop(indx)
  #peak1_val = 

#for indx in idx_dict["Index"]:
 # if indx not in 
  #print(indx)
#hash_map_pool_1["Concentration".pop()]

#print(hash_map_pool_1)
#for idx_list in idx_dict["Index"]:
  #for xyz in idx_list:
   # print(xyz)
    #hash_map_pool_1["Compiled_Concentration"].append(hash_map_pool_1["Concentration"][xyz])
    

    #print(hash_map_pool_1["Compiled_Concentration"])    


    
    #pool_1["Size"].append()
   
    #find_peak = range)


          #### 






  

    
  



# remove all values in this occurence within all hash_map_pool keys

# pool 1
#for key in hash_map_pool_1.keys(): # key = size, concentration, marker
    #for x in hash_map_pool_1[key]:
      #print('NaN')
    #print(x)

# peaks = []  # values with just peaks
# for num, x in enumerate(hash_map_pool_1["Markers"]):
#     if num not in idx_remove: # pulling out desired peaks
#         peaks.append(num)

#print(peaks)   # all the lines with peaks



#print(os.getcwd())
#
#Molarity = moles of solute/volume = mol/m^3
  



#print(hash_map_pool_1.keys()) 



# defined list: hash_map_pool_1["Markers"]

# find the index with Lower Marker
#idx_marker = hash_map_pool_1["Size"].index()
#print(idx_marker)

#replace lower marker value with ""
#hash_map_pool_1=hash_map_pool_1[:"Markers"]+[]
#for tracker, val in enumerate(hash_map_pool_1["Markers"]): 
  #print('x')
  


      
# Next Step: P2/P1 Ratio
