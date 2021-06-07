#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 19:40:36 2020

@author: elizabethgonzalez
"""
import numpy as np
import pandas as pd

G = 6.6743 * 10**-11
solar_mass = 2 * 10**30
earth_mass = 6 *10**24
AU = 1.5 *10**11
solar_radius = 7 * 10**8

def metadata(filename):
    for i in range(36):
        f.readline()
    
#use their timestamp & error
t = open("hip5364.vels", 'r')
lines=t.readlines()
time=[]
error = []
for x in lines:
    time.append(float(x.split('    ')[1]))
    error.append(float(x.split('    ')[3]))
t.close()
time = np.array(time)

        
#generating sine waves
f = open("ajaaf477t2_mrt copy.txt", 'r')
metadata(f)


#read in table
data = {}
for elem in range(5391):
    line = f.readline().strip().split()
    if line[0] not in data:
        data[str(line[0])] = []
    data[str(line[0])].append(line[0:])
    
#find individual planet RV sine wave
phase_list = [[]for i in range(953)]
true_amp = []
true_periods = []
transit = []
detect = []
ID = []
planet = []
total_planet = []
planet_mass = []
major_axis = []
for i in range(1,953):
    x = str(i)
    tot = len(data[x])
    tot_planet = [len(data[x])]*tot
    total_planet.extend(tot_planet)
    for j in range(0,len(data[x])):
        e = float(data[x][j][-4])
        true_periods.append(float(data[x][j][8]))
        transit.append(float(data[x][j][-2]))
        detect.append(float(data[x][j][-1]))
        ID.append(float(data[x][j][0]))
        planet.append(float(data[x][j][1]))
        planet_mass.append(float(data[x][j][11]))
        period = float(data[x][j][8])*86400
        period_sin = float(data[x][j][8])
        inc = float(data[x][j][-8]) * (np.pi/180)
        Teff = int(data[x][j][2])
        p_mass = float(data[x][j][-7])
        major_axis.append(float(data[x][j][-9]))
        
        #finding star mass & radius
        star_mass = -22.296508 + (1.5446387*10**-2)*Teff - (3.488452*10**-6)*(Teff**2) + (2.64961*10**-10)*(Teff**3)
        star_radius = -16.883175 + 1.1835396*10**-2*Teff - (2.70872*10**-6)*(Teff**2) + 2.105028*10**-10*(Teff**3)
     
        p_mass = p_mass * earth_mass
        star_mass = star_mass * solar_mass
        star_radius = star_radius * solar_radius       
        
        #finding semi-major axis both ways
        a_r = float(data[x][j][-9])
        a_table = star_radius*a_r
        
        #generating sine wave
        phase = np.random.uniform(0,2*np.pi)
        amp_table = G**(1/2)*(star_mass + p_mass)**(-1/2)*(a_table)**(-1/2)*p_mass*np.cos(inc)
        RV = amp_table * np.sin(((time*2*np.pi)/period_sin - phase))
        amp = G**(1/2)*(star_mass + p_mass)**(-1/2)*(a_table)**(-1/2)*p_mass*np.cos(inc)
        
        true_amp.append(amp)
        data[x][j] = RV
        phase_list[i].append(phase)

#add individual planet sine wave to get total star RV
rv_sum = [[]for i in range(953)] #lst[0] will be blank therefore index = star id 
for i in range(1, 953):
    x = str(i)
    rv_sum[i] = np.add(data[x][0], data[x][1])
    if len(data[x]) > 2:
        for j in range(2, len(data[x])):
            rv_sum[i] = np.add(rv_sum[i], data[x][j])
    rv_sum[i] = np.add(rv_sum[i], error)   
  

for i in range(1, 953):      
    rv_sum[i] = [round(num, 1) for num in rv_sum[i]]
    with open("RV_Star{}.txt".format(i), 'w') as f:
        for a,b,c in zip(time, rv_sum[i], error):
            f.write('{0:15}{1:15}{2:15}\n'.format(a, b, c)) 
    

title = ['Star', 'Total_Planet', 'True_Period', 'True_Amp', 'Transit', 'Detect', 'Planet_Mass', 'Major_Axis']     

df = pd.DataFrame(columns=title)
df['Star'] = ID
df['Total_Planet'] = total_planet
df['True_Period'] = true_periods
df['True_Amp'] = true_amp
df['Transit'] = transit
df['Detect'] = detect
df['Planet_Mass'] = planet_mass
df['Major_Axis'] = major_axis

df.to_pickle("true_data.pkl")