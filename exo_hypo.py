#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 19:12:57 2020

@author: elizabethgonzalez
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
df = pd.read_pickle("found_data.pkl")
df_true = pd.read_pickle("true_data.pkl")
#print(df, "\n", df_true)

df_true["Exo_Detection"] = np.nan
df_true["Found_Period"] = np.nan
df_true["Period_Error"] = np.nan
df_true["Found_Amp"] = np.nan
df_true["Amp_Error"] = np.nan


for index, row in df_true.iterrows():
    star = df_true.iloc[index]['Star']
    period = df_true.iloc[index]['True_Period']
    exostriker = df.loc[df['Star'] == star]
    for i in range(0, len(exostriker)):
        found_per = exostriker.iloc[i]['Found_Period']
        if found_per < period*1.05 and found_per > period*0.95:
            df_true.at[index, "Exo_Detection"] = 1 #found
            df_true.at[index, "Found_Period"] = exostriker.iloc[i]['Found_Period']
            df_true.at[index, "Period_Error"] = exostriker.iloc[i]['Period_Error']
            df_true.at[index, "Found_Amp"] = exostriker.iloc[i]['Found_Amp']
            df_true.at[index, "Amp_Error"] = exostriker.iloc[i]['Amp_Error']

pd.set_option("display.max_rows", None, "display.max_columns", None)
#print(df_true)

for i in range(1,953):
    x = df_true[df_true['Star'] == i]
    for j in range(0,len(x)):
        if x.iloc[j]['Exo_Detection'] == 1.0: #blue for detected
            plt.plot(np.log(x.iloc[j]['Major_Axis']),x.iloc[j]['Star'],'.', color = 'b', markersize=x.iloc[j]['Planet_Mass'])
        if x.iloc[j]['Exo_Detection'] != 1: #red for notdetected
            plt.plot(np.log(x.iloc[j]['Major_Axis']),x.iloc[j]['Star'],'.', color = 'r', markersize=x.iloc[j]['Planet_Mass'])
plt.show()
          
'''
#determining if neighbor planets are detected for those that wrongly found
number_neighbors = []
number_undetected_neighbors = []
one_undetected_neighbor = []
two_undetected_neighbor = []
two_detected_neighbor = []
one_detected_neighbor = []
one_undetected_ofone = []
for i in range(1, 953):
    x = df_true[df_true['Star'] == i]
    for j in range(0, len(x)):
        if x.iloc[j]['Exo_Detection'] != 1.0:
            if j == 0:
                if df_true.iloc[j + 1]['Exo_Detection'] != 1:
                    one_undetected_ofone.append(1) 
                    number_neighbors.append(1)
                    number_undetected_neighbors.append(1)
                if df_true.iloc[j + 1]['Exo_Detection'] == 1:
                    one_detected_neighbor.append(1)  
                    number_neighbors.append(1)
                    number_undetected_neighbors.append(0)
            if j == (len(x)-1):
                if df_true.iloc[j - 1]['Exo_Detection'] != 1:
                    one_undetected_ofone.append(1)
                    number_neighbors.append(1)
                    number_undetected_neighbors.append(1)
                if df_true.iloc[j - 1]['Exo_Detection'] == 1:
                    one_detected_neighbor.append(1)
                    number_neighbors.append(1)
                    number_undetected_neighbors.append(0)
            if j != 0 and j != (len(x)-1):
                if df_true.iloc[j - 1]['Exo_Detection'] != 1 and df_true.iloc[j + 1]['Exo_Detection'] != 1:
                    two_undetected_neighbor.append(1)
                    number_undetected_neighbors.append(2)
                    number_neighbors.append(2)
                if df_true.iloc[j - 1]['Exo_Detection'] == 1 and df_true.iloc[j + 1]['Exo_Detection'] == 1:
                    two_detected_neighbor.append(1)
                    number_neighbors.append(2)
                    number_undetected_neighbors.append(0)
                if df_true.iloc[j - 1]['Exo_Detection'] != 1 and df_true.iloc[j + 1]['Exo_Detection'] == 1:
                    one_undetected_neighbor.append(1)
                    number_neighbors.append(2)
                    number_undetected_neighbors.append(1)
                if df_true.iloc[j - 1]['Exo_Detection'] == 1 and df_true.iloc[j + 1]['Exo_Detection'] != 1:
                    one_undetected_neighbor.append(1) 
                    number_neighbors.append(2)
                    number_undetected_neighbors.append(1)
        if x.iloc[j]['Exo_Detection'] == 1.0:
            if j == 0:
                if df_true.iloc[j + 1]['Exo_Detection'] != 1:
                    number_neighbors.append(1)
                    number_undetected_neighbors.append(1)
                if df_true.iloc[j + 1]['Exo_Detection'] == 1: 
                    number_neighbors.append(1)
                    number_undetected_neighbors.append(0)
            if j == (len(x)-1):
                if df_true.iloc[j - 1]['Exo_Detection'] != 1:
                    number_neighbors.append(1)
                    number_undetected_neighbors.append(1)
                if df_true.iloc[j - 1]['Exo_Detection'] == 1:
                    number_neighbors.append(1)
                    number_undetected_neighbors.append(0)
            if j != 0 and j != (len(x)-1):
                if df_true.iloc[j - 1]['Exo_Detection'] != 1 and df_true.iloc[j + 1]['Exo_Detection'] != 1:
                    number_undetected_neighbors.append(2)
                    number_neighbors.append(2)
                if df_true.iloc[j - 1]['Exo_Detection'] == 1 and df_true.iloc[j + 1]['Exo_Detection'] == 1:
                    number_neighbors.append(2)
                    number_undetected_neighbors.append(0)
                if df_true.iloc[j - 1]['Exo_Detection'] != 1 and df_true.iloc[j + 1]['Exo_Detection'] == 1:
                    number_neighbors.append(2)
                    number_undetected_neighbors.append(1)
                if df_true.iloc[j - 1]['Exo_Detection'] == 1 and df_true.iloc[j + 1]['Exo_Detection'] != 1:
                    number_neighbors.append(2)
                    number_undetected_neighbors.append(1)
df_true['Number_Neighbors'] = number_neighbors
df_true['Number_UndetectedNeighbors'] = number_undetected_neighbors

#number found out of total
detected = df_true.loc[df_true['Exo_Detection'] == 1.0]
undetected = df_true.loc[df_true['Exo_Detection'] != 1.0]
print("{} found out of {} are detected".format(len(detected), len(df_true)))
print('\n')

print("For planets with two neighbors: ")
print("{} out of {} found have 1 undetected/detected neighbor".format(len(one_undetected_neighbor), len(undetected)))
print("{} out of {} found have 2 undetected neighbors".format(len(two_undetected_neighbor), len(undetected)))
print("{} out of {} found have 2 detected neighbors".format(len(two_detected_neighbor), len(undetected)))            
print("For planets with one neighbor: ")
print("{} out of {} found have 1 detected neighbors".format(len(one_detected_neighbor), len(undetected))) 
print("{} out of {} found have 1 undetected neighbors".format(len(one_undetected_ofone), len(undetected)))

detected_withall_neighbors = []
detected_withsome_neighbors = []
for i in range(0,len(detected)):
    if detected.iloc[i]['Number_UndetectedNeighbors'] == 0:
        detected_withall_neighbors.append(1)
    else:
        detected_withsome_neighbors.append(1)
print("{} out of {} detected planets have all neighbors detected".format(len(detected_withall_neighbors), len(detected)))        
print("{} out of {} detected planets have 1+ neighbors undetected".format(len(detected_withsome_neighbors), len(detected)))
      
undetected_withall_neighbors = []
undetected_withsome_neighbors = []
for i in range(0,len(undetected)):
    if undetected.iloc[i]['Number_UndetectedNeighbors'] == 0:
        undetected_withall_neighbors.append(1)
    else:
        undetected_withsome_neighbors.append(1)
print("{} out of {} undetected planets have all neighbors detected".format(len(undetected_withall_neighbors), len(undetected)))        
print("{} out of {} undetected planets have 1+ neighbors undetected".format(len(undetected_withsome_neighbors), len(undetected)))    

plt.scatter(df_true['Planet_Mass'], (df_true['Number_Neighbors'] - df_true['Number_UndetectedNeighbors']))
plt.title("Detected Neighbor Ratio vs Planet Mass")
plt.xlabel("Mass")
plt.ylabel("Detected Neighbor Ratio")
plt.show() 

#found amplitude error between true and found amp
error_amp = []
for i in range(0,len(df_true)):
    true = float(df_true.iloc[i]['True_Amp'])
    found = float(df_true.iloc[i]['Found_Amp'])
    err = (true-found)/true
    error_amp.append(err)
plt.hist(error_amp,100)
plt.title("Found Amplitude Error")
plt.xlabel("Error")
plt.ylabel("Frequency")
plt.show()
df_true['error_amp'] = error_amp

#found period error between true and found period
error_period = []
for i in range(0,len(df_true)):
    true = float(df_true.iloc[i]['True_Period'])
    found = float(df_true.iloc[i]['Found_Period'])
    err = (true-found)/true
    error_period.append(err)
plt.hist(error_period,100)
plt.title("Found Period Error")
plt.xlabel("Error")
plt.ylabel("Frequency")
plt.show()
df_true['error_period'] = error_period

detected = df_true.loc[df_true['Exo_Detection'] == 1.0]
undetected = df_true.loc[df_true['Exo_Detection'] != 1.0]
print("Median Amplitude: {}".format(np.median(detected['error_amp'])))
print("Median Period: {}".format(np.median(detected['error_period'])))

#plot period vs amp and color code if found or not
detected = df_true.loc[df_true['Exo_Detection'] == 1.0]
plt.scatter(np.log10(detected['True_Period']), detected['True_Amp'], color='k', label = 'Detected')
plt.scatter(np.log10(undetected['True_Period']), undetected['True_Amp'], color='g', label = 'Undetected')
plt.title("True Period v True Amplitude")
plt.xlabel("Period")
plt.ylabel("Amplitude")
plt.legend()
plt.show() 
'''


#df_true.to_pickle("complete_table.pkl")   