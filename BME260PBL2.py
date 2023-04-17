import numpy as np
import matplotlib.pyplot as plt

#%% Healthy Steady State 

#%% Known Values for Steady State:
    
n1 = 266 / 86400 #nmol/day / 864000 seconds/day = nmol/sec
vitKa = 1.45 #nmol/L 
Xa = 170 #nmol/L 
proth = 1400 #nmol/L
fbng = 9000 #nmol/L
thr = 1 #nmol/L
tf1 = 1.6*10**-5 #nmol/L (healthy)
X = 160 #nmol/L
Xi = 160 #nmol/L
vitKi = vitKa
fbr = fbng
Vol_blood = 5 #Liters
Vol_liver = Vol_blood*0.125

#rate constants
kcm = 1559 #nmol/L
kc2 = 55.9 #1/s
kd2m = 1600 #nmol/L
kd22 = 57 #1/s
kf1 = 0.025 #L/(s*nmol)
kem = 7200 #nmol/L
ke2 = 84 #1/s


#%% Box D (Thrombin Activation)

n8 = kd22*Xa*proth/(proth+kd2m) * Vol_liver #prothrombin
n11 = kd22*Xa*proth/(proth+kd2m) * Vol_liver #thrombin

#%% Box E (Capillaries)

n9 = ((ke2*thr*fbng)/(fbng+kem)) * Vol_blood #fibrinogen
n16 = ((ke2*thr*fbng)/(fbng+kem)) * Vol_blood #fibrin
k_thr_deg = n11/(thr*Vol_blood) #thrombin #1/s
n10 = k_thr_deg * thr * Vol_blood #nmol/(L * s)

#%% Box A (Intestine)

n2 = n1*0.45
n3 = n1*0.55
n1 = n3 + n2

#%% Box F (Injury Site)

n12 = kf1*tf1*X*Vol_blood #X #nmol/(L*s)
n13 = kf1*tf1*X*Vol_blood #Xi #nmol/(L*s)
k_fbr_deg = n16/(fbr*Vol_blood) #fibrin #1/s
n14 = k_fbr_deg * fbr * Vol_blood

#%% Box D (Intermediate Equation)
#solving for rate constant for this

kd1 = (n13/Vol_liver)/(vitKa*Xi) #Xi #L/(nmol*s)
n7 = kd1*vitKa*Xi *Vol_liver#vitKi #nmol/(L*s)
n6 = kd1*vitKa*Xi *Vol_liver#vitKa #nmol/(L*s)
k_Xa_deg = (kd1*vitKa*Xi) / Xa #Xa #1/s
n15 = k_Xa_deg * Xa * Vol_liver

#%% Box B (Liver Production Box)

n4 = 0.55*n1
fbng_gen = n9  #fibrinogen
proth_gen = n8 #prothrombin
X_gen = n12 #X
n5 = n3 - n4 + n7 #vitKi

#%% Box C (Vitamin K Reduction)

VKORC1 = (n5/Vol_liver)*((vitKi + kcm)/((vitKi) * kc2))


print("n1: {:.3e}".format(n1))
print("n2: {:.3e}".format(n2))
print("n3: {:.3e}".format(n3))
print("n4: {:.3e}".format(n4))
print("n5: {:.3e}".format(n5))
print("n6: {:.3e}".format(n6))
print("n7: {:.3e}".format(n7))
print("n8: {:.3e}".format(n8))
print("n9: {:.3e}".format(n9))
print("Thrombin Degredation, n10: {:.3e}".format(n10))
print("n11: {:.3e}".format(n11))
print("n12: {:.3e}".format(n12))
print("n13: {:.3e}".format(n13))
print("n14: {:.3e}".format(n14))
print("n15: {:.3e}".format(n15))
print("n16: {:.3e}".format(n16))
print("fbng_gen: {:.3e}".format(fbng_gen))
print("proth_gen: {:.3e}".format(proth_gen))
print("X_gen: {:.3e}".format(X_gen))
print("kd1: {:.3e}".format(kd1))
print("VKORC1: {:.3e}".format(VKORC1))
print("K_thrombin_deg: {:.3e}".format(k_thr_deg))
print("K_firbin_deg: {:.3e}".format(k_fbr_deg))
print("K_Xa_deg: {:.3e}".format(k_Xa_deg))


#%%Graphing Data

fig1 = plt.figure(num=1, clear=True)
ax1 = fig1.subplots(1,1)
compound = ['Fibrinogen','Fibrin', 'Thrombin','Prothrombin']
flow_rate = [n9,n16, n11, n8]
plt.bar(compound, flow_rate)
