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
sig = 1

#rate constants
dk = 57                 #1/s
fk = 0.025              #L/(s*nmol)
ekm = 7200              #nmol/L
ek2 = 84                #1/s

#%% Box D (Thrombin Activation)

n6 = dk*Xa*vitKa *Vol_liver
n12 = dk*Xa*vitKa *Vol_blood
n7 = dk*Xa*vitKa *Vol_liver
n10 = dk*Xa*vitKa *Vol_liver#thrombin

#%% Box F (Injury site)

R = n12/sig
n11 = sig*R *Vol_liver
ngen = sig*R #TF generation

#%%Box E (Capillaries)

n8 = (ek2*thr*fbng)/(fbng + ekm) * Vol_blood  #fibrinogen
n15 = (ek2*thr*fbng)/(fbng + ekm) * Vol_blood #fibrin
ek = n10/thr
n9 = n11 *Vol_blood

#Box F part 2
fk = n15/fbr
n13 = fk*fbr *Vol_blood


#%% Box A (Intestine)

n2 = n1*0.45
n3 = n1*0.55
n1 = n3 + n2


#%% Box B (Liver Production)

n4 = 0.55*n1
fbng_gen = n8  #fibrinogen
X_gen = n12 #X
n5 = n3 - n4 + n7 #vitKi

#%% Box C (Vitmain K reduction)

ck = n5/vitKi #vitki
n6 = ck*vitKi *Vol_blood#vitKa

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
print("n14: {:.3e}".format(n15))
print("fbng_gen: {:.3e}".format(fbng_gen))
print("X_gen: {:.3e}".format(X_gen))
print("kd1: {:.3e}".format(dk))
print("K_thrombin_deg: {:.3e}".format(ek))
print("K_firbin_deg: {:.3e}".format(fk))

#%%Graphing Data

#fig1 = plt.figure(num=1, clear=True)
#ax1 = fig1.subplots(1,1)
#compound = ['Fibrinogen','Fibrin', 'Thrombin','Prothrombin']
#flow_rate = [n9,n16, n11]
#plt.bar(compound, flow_rate)