import numpy as np
import matplotlib.pyplot as plt

#%% Healthy Steady State 

#%% Known Values for Steady State:
n1 = 266 / 86400        #nmol/day / 864000 seconds/day = nmol/sec
vitKa = 1.45            #nmol/L 
Xa = 170                #nmol/L 
fbng = 9000             #nmol/L
thr = 1                 #nmol/L
tf = 1.6*10**-5        #nmol/L (healthy)
X = 160                 #nmol/L
vitKi = vitKa
fbr = fbng
Vol_blood = 5           #Liters
Vol_liver = Vol_blood*0.125

#rate constants
kcm = 1559              #nmol/L
kc2 = 55.9              #1/s
kd2m = 1600             #nmol/L
kdd2 = 57               #1/s
kf1 = 0.025             #L/(s*nmol)
kem = 7200              #nmol/L
ke2 = 84                #1/s

#%% Box D (Thrombin Activation)

n11 = kdd2*Xa/() * Vol_liver #thrombin

#%%Box E (Capillaries)

n6 = (ek2*thr*fbng)/(fbng + ekm) #fibrinogen
n16 = (ek2*thr*fbng)/(fbng + ekm) #fibrin

#%% Box A (Intestine)


#%% Box F (Injury site)


#%% Box B (Liver Production)

#%% Box C (Vitmain K reduction)


#%%Graphing Data

fig1 = plt.figure(num=1, clear=True)
ax1 = fig1.subplots(1,1)
compound = ['Fibrinogen','Fibrin', 'Thrombin','Prothrombin']
flow_rate = [n9,n16, n11, n8]
plt.bar(compound, flow_rate)