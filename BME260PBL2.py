#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  9 19:59:14 2023

@author: jameskoconis
"""

import numpy as np
import matplotlib.pyplot as plt

#%% Healthy Steady State 

#%% Known Values for Steady State:
    
n1 = 266 #nmol/day
vitKa = 1.45 #nmol/L
Xa = 170 #nmol/L
proth = 1400 #nmol/L
fbng = 8823 #nmol/L
thr = 1 #nmol/L
tf1 = 1.6*10**5 #nmol/L (healthy)
X = 160 #nmol/L
Xi = 160 #nmol/L
vitKi = vitKa

#rate constants
kcm = 1950 #nmol/L
kd2m = 1600 #nmol/L
kd22 = 57 #1/s
kf1 = 0.025 #L/(s*nmol)
kem = 7200 #nmol/L
ke2 = 84 #1/s


#%% Box D (Thrombin Activation)

n8 = kd22*Xa*proth/(proth+kd2m) #prothrombin
n11 = kd22*Xa*proth/(proth+kd2m) #thrombin

#%% Box E (Capillaries)

n9 = ke2*thr*fbng/(fbng+kem) #fibrinogen
n16 = ke2*thr*fbng/(fbng+kem) #fibrin
n10 = n11 #thrombin

#%% Box A (Intestine)

n2 = n1*0.45
n3 = n1*0.55
n1 = n3 + n2

#%% Box F (Injury Site)

n12 = kf1*tf1*X #X #nmol/(L*s)
n13 = kf1*tf1*X #Xi #nmol/(L*s)
n14 = n16 #fibrin #nmol/(L*s)

#%% Intermediate Equation (part of Box D)
#solving for rate constant for this

kd1 = n13/(vitKa*Xi) #Xi #L/(nmol*s)
n7 = kd1*vitKa*Xi #vitKi #nmol/(L*s)
n6 = kd1*vitKa*Xi #vitKa #nmol/(L*s)
n15 = kd1*vitKa*Xi #Xa #nmol/(L*s)
n13 = kd1*vitKa*Xi #Xi #nmol/(L*s)

#%% Liver Production Box

n4 = 0.55*n1
fbng_gen = n9  #fibrinogen
proth_gen = n8 #prothrombin
X_gen = n10 #X
n5 = n3 - n4 + n7 #vitKi

#%% Box C (Vitamin K Reduction)

kc2_VKORC1 = n5*(vitKi + kcm)/vitKi


print("n1", n1)
print("n2", n2)
print("n3", n3)
print("n4", n4)
print("n5", n5)
print("n6", n6)
print("n7", n7)
print("n8", n8)
print("n9", n9)
print("n10", n10)
print("n11", n11)
print("n12", n12)
print("n13", n13)
print("n14", n14)
print("n15", n15)
print("fbng_gen", fbng_gen)
print("proth_gen", proth_gen)
print("X_gen", X_gen)
print("kd1", kd1)
print("kc2_VKORC1", kc2_VKORC1)


#%%Graphing Data

fig1 = plt.figure(num=1, clear=True)
ax1 = fig1.subplots(1,1)
compound = ['Fibrinogen','Fibrin', 'Thrombin','Prothrombin']
flow_rate = [n9,n16, n11, n8]
plt.bar(compound, flow_rate)

#%% Healthy Dynamic

#%% Known Values

#rate constants
kcm = 1950 #nmol/L
kd2m = 1600 #nmol/L
kd22 = 57 #1/s
kf1 = 0.025 #L/(s*nmol)
kem = 7200 #nmol/L
ke2 = 84 #1/s
kd1 = 2758.62 #L/(nmol*s)
kc2_VKORC1 = 861329655 #(nmol/L)**2 *1/s

n1 = 266 #nmol/day
vitKa = 1.45 #nmol/L
Xa = 170 #nmol/L
proth = 1400 #nmol/L
fbng = 8823 #nmol/L
thr = 1 #nmol/L
tf1 = 1.6*10**5 #nmol/L (healthy)
X = 160 #nmol/L
Xi = 160 #nmol/L
vitKi = vitKa