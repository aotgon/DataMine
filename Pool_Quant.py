# Goal: Determine for potential overloading of cells due to too many PCR cycles

# initialize packages 
from asyncore import write
from multiprocessing import pool
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
  pool_1["Concentration_Peak"].append(pks)
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
  pool_2["Concentration_Peak"].append(pks)
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
  pool_3["Concentration_Peak"].append(pks)
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
  pool_4["Concentration_Peak"].append(pks)
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
fieldnames_1 = ['Pool1_Concentration_ng/uL', 'Pool1_Size_bp',
'Pool2_Concentration_ng/uL ', 'Pool2_Size_bp',
'Pool3_Concentration_ng/uL ', 'Pool3_Size_bp',
'Pool4_Concentration_ng/uL ', 'Pool4_Size_bp',
]
#sub_fieldname = ['Pool_Concentration','Pool_Size']


# find how many peaks in each pool
#peak_qty = len(str(pool_4["Size"][0][0]))
 # construct dataframe from a dictionary

#peak_qty = len(str(pool_1["Size"][0][0]))
temp_list = []
df = []
for x, y in enumerate(pool_1["Size"]): # number of peaks
  #print(x,y) # x is idx, y is the size 
  #idk =  (pool_1["Concentration_Peak"][x]) #[list:{x,y}]
  #idk_1 = idk[0] # [x]
  #print(idk_1)

  current_pk_conc = pool_1["Concentration_Peak"][x] #[1.58],[1.54],[1.56]
  concen_only = current_pk_conc[0] # only first component of list
  #print(concen_only)

  current_pk_size = pool_1["Size"][x]
  size_only = current_pk_size[0]
  
  #print(size_only)
  #print(pool_1["Size"][x][0])
  d = {'Concentration_Peak': concen_only,
      'Size':size_only,
       }

  #temp_df = pd.DataFrame(data=d,ignore_index = True)
  
  #df.append(temp_df)
  #temp_df = pd.DataFrame()
  
df = pd.DataFrame(data=pool_1)
#df.set_index('Concentration_Peak', 1**peak_qty)
print(df)
  

#df = pd.DataFrame(data=pool_1)
#print(df)

######################################################################
# #subheader = ['Concentration_Peak_ng/ul','','Size_bp']
# with open('/Users/aotgonbayar/Desktop/Python/code.csv','w') as f:# open the file in the write mode

# # # create the csv writer

#   writer_csv = csv.DictWriter(f, fieldnames= fieldnames_1)

# # create subheader
#   writer_csv.writeheader()   # dataframe, 

#   for sub_fieldname in fieldnames_1:
    
#     if 'Pool1_Concentration' in sub_fieldname:
#       #counter= 0
#       for idx_current, key_current in enumerate(pool_1["Concentration_Peak"]):
#         #counter+=1
#         #writer_csv.writerow({'Pool1_Concentration_ng/uL':pool_1["Concentration_Peak"[idx_current][0]]})
#         #print(pool_1["Concentration_Peak"][idx_current][0])
#         writer_csv.writerow({'Pool1_Concentration_ng/uL':pool_1["Concentration_Peak"][idx_current][0],
#                             'Pool1_Size_bp':pool_1["Size"][idx_current][0]})
#     elif 'Pool1_Concentration' in sub_fieldname:
#       for idx_current, key_current in enumerate(pool_2["Concentration_Peak"]):
#         #counter+=1
#         #writer_csv.writerow({'Pool1_Concentration_ng/uL':pool_1["Concentration_Peak"[idx_current][0]]})
#         #print(pool_1["Concentration_Peak"][idx_current][0])
#         writer_csv.writerow({'Pool2_Concentration_ng/uL':pool_2["Concentration_Peak"][idx_current][0],
#                             'Pool2_Size_bp':pool_2["Size"][idx_current][0]})                      
#############################################################################
        #writer_csv.writerow({'Pool1_Size_bp':pool_1["Size"][idx_current][1]})

        #print(idx_current,key_current) [0,1,2],[1.58,0.0332]
      
      #print(sub_fieldname)
      #writer_csv.writerow(pool_1["Concentration_Peak"])
        #print(idx_current,key_current)
      #writer_csv.writerow(sub_fieldname,{pool_1["Concentration_Peak"][[0][0]]})
        #writer_csv.writerow({pool_1["Concentration_Peak"][idx_current][0]})





  
  



  #for x_idx, key_x in enumerate(pool_1["Concentration_Peak"]):
    #writer_csv.writerow({'Pool1':pool_1["Concentration_Peak"][x_idx][0]})

  

# # write a row to the csv file
  #writer_csv.writerow(header)
  
  #writer_csv.writerow(pool_1["Concentration_Peak"])

  #for sub_head in header:
    #if 'Pool' in sub_head: # if header has pool in their name: then make a subheader
      #writer_csv.writerow(subheader) # subheader with concentration and size (1x2 column)
    #for key in pool_1["Concentration_Peak"]: 
      #writer_csv.writerow([key[0],key[1]])

  #writer_csv.writerow([pool_1[],pool_2,pool_3,pool_4])






  

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
