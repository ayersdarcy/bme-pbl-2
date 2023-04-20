import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
mpl.rcParams['font.size'] = 14

# %% Initial Values
#reaction for vitamin k reduction (box c)
#rate constants
ck = 1.7*10
#reaction for thrombin activation (box d)
#rate constants
kd1 = 0.1
#reaction for injury site (box f)
#rate constants
fk = 0.025
#reaction for fibrin gen (box e)
#rate constants
ek1 = 0.1
ek_1 = 63.6
ek2 = 84
#degradation rate constants
k_thrombin_deg = 3.081*10**-2    # 1/sec
k_fibrin_deg = 5.185*10**-4     # 1/sec
k_vitKa_exc = 1.868*10**-3

#Reaction Rate
R = 1.232*10**2

#Streams
n1=  0.0030787037   #nmol/sec
n2 = 1.385*10**-3   #nmol/sec
n3 = 1.693*10**-3   #nmol/sec
n4 = 1.693*10**-3   #nmol/sec
n5 = 1.541*10     #nmol/sec
n6 = 1.541*10     #nmol/sec
n7 = 1.541*10     #nmol/sec
n8 = 2.333*10**2    #nmol/sec
n9 = 1.541*10    #nmol/sec
n10 = 1.541*10   #nmol/sec
n11 = 7.703*10   #nmol/sec
n12 = 1.232*10**2    #nmol/sec
n13 = 2.333*10**2    #nmol/sec
n14 = 2.333*10**2   #nmol/sec

#known values for steady state
vitKa = 1.45                #nmol/L
Xa = 170                    #nmol/L
fbng = 9000                 #nmol/L
thr = 1                    #nmol/L
tf1 = 1.6*10**-4            #nmol/L
X = 160                     #nmol/L
vitKi = vitKa                #nmol/L
fbr = fbng                  #nmol/L
Vol_blood = 5               #Liters
Vol_liver = Vol_blood*0.125 #Liters
ES = 0


#rate constants
kc2 = 3*10**6 #1/sec
kd22 = 57 #1/sec
kf1 = 0.025 #L/sec*nnmol
ke2= 84 #1/sec

# %% Calculations
tspan = np.linspace(0, 60, 10000) #time point for each second 
y0 = np.array([n5, n6, n7, n10, n12, n14,
                vitKa, vitKi, X, Xa, thr, 
                fbng, ES, tf1, fbr])

def odefunc(y0, t): 
    n5 = y0[0]
    n6 = y0[1]
    n7 = y0[2]
    n10 = y0[3]
    n12 = y0[4]
    n14 = y0[5]
    vitKa = y0[6]
    vitKi = y0[7]
    X = y0[8]
    Xa = y0[9]
    thr = y0[10]
    fbng = y0[11]
    ES = y0[12]
    tf1 = y0[13]
    fbr = y0[14]

    #box D
    dVitKa = (n6 * Vol_liver) - (kd1 * vitKa * Xa)
    dn7 = kd1 * vitKa * Xa * Vol_liver
    dn10 = dn7
    dXa = n12 / Vol_liver - kd1*vitKa*Xa

    #box E
    dthr = (-ek1*thr*fbng) + (ek_1*ES) + (ek2*ES) - (k_thrombin_deg*thr) + (n10/Vol_blood)
    dfbng = (-ek1*thr*fbng) + (ek_1*ES) + (n8/Vol_blood)
    dES = (ek1*thr*fbng) - (ek_1*ES) - (ek2*ES)
    dn14 = (ek2*ES*Vol_blood)

    #box B
    dVitKi = (n3/Vol_liver) - (k_vitKa_exc*vitKi) - (n5/Vol_liver) + (n7/Vol_liver)

    #box C
    dn5 = (ck*vitKi*Vol_liver)
    dn6 = (ck*vitKi*Vol_liver)

    #box F
    #dTF = (-5/10**6)*R
    dTF = -tf1*0.1
    dX = (n11/Vol_blood) - (R/Vol_blood)
    dn12 = R
    dfbr = (n14/Vol_blood) - k_fibrin_deg*fbr

    return np.array([dn5-n5, dn6-n6, dn7-n7, dn10-n10, dn12-n12, dn14-n14,
                     dVitKa, dVitKi, dX, dXa, dthr, 
                     dfbng, dES, dTF, dfbr])

output = odeint(odefunc, y0, tspan)

vitKaResult = output[:,6]
XResult = output[:,8]
XaResult = output[:,9]
thrResult = output[:,10]
fbrResult = output[:,14]
n6Result = output[:,1]


# CHANGE TO DISEASED DYNAMIC 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp

# %% Initial Values
#reaction for vitamin k reduction (box c)
#rate constants
ck = 1.7*10
#reaction for thrombin activation (box d)
#rate constants
kd1 = 0.1
#reaction for injury site (box f)
#rate constants
fk = 0.025
#reaction for fibrin gen (box e)
#rate constants
ek1 = 0.1
ek_1 = 63.6
ek2 = 84
#degradation rate constants
k_thrombin_deg = 3.081*10**-2    # 1/sec
k_fibrin_deg = 5.185*10**-4     # 1/sec
k_vitKa_exc = 1.868*10**-3

#Reaction Rate
R = 1.232*10**2

#Streams
n1=  1.539*10**-3   #nmol/sec
n2 = 6.927*10**-4   #nmol/sec
n3 = 8.466*10**-4   #nmol/sec
n4 = 8.466*10**-4   #nmol/sec
n5 = 7.703          #nmol/sec
n6 = 7.703          #nmol/sec
n7 = 7.703          #nmol/sec
n8 = 2.333*10**2    #nmol/sec
n9 = 7.703          #nmol/sec
n10 = 7.703         #nmol/sec
n11 = 3.852*10   #nmol/sec
n12 = 6.162*10    #nmol/sec
n13 = 2.333*10**2    #nmol/sec
n14 = 2.333*10**2   #nmol/sec

#known values for steady state
vitKa = 0.725               #nmol/L
Xa = 170                    #nmol/L
fbng = 9000                 #nmol/L
thr = 1                    #nmol/L
tf1 = 1.6*10**-4            #nmol/L
X = 160                     #nmol/L
vitKi = vitKa                #nmol/L
fbr = fbng                  #nmol/L
Vol_blood = 5               #Liters
Vol_liver = Vol_blood*0.125 #Liters
ES = 0


#rate constants
kc2 = 3*10**6 #1/sec
kd22 = 57 #1/sec
kf1 = 0.025 #L/sec*nnmol
ke2= 84 #1/sec

# %% Calculations
y0 = np.array([n5, n6, n7, n10, n12, n14,
                vitKa, vitKi, X, Xa, thr, 
                fbng, ES, tf1, fbr])

def odefunc2(y0, t): 
    #initialize all the values that change at each time point 
    #some of the stream values 
    n5 = y0[0]
    n6 = y0[1]
    n7 = y0[2]
    n10 = y0[3]
    n12 = y0[4]
    n14 = y0[5]
    vitKa = y0[6]
    vitKi = y0[7]
    X = y0[8]
    Xa = y0[9]
    thr = y0[10]
    fbng = y0[11]
    ES = y0[12]
    tf1 = y0[13]
    fbr = y0[14]

    #box D
    dVitKa = (n6 * Vol_liver) - (kd1 * vitKa * Xa)
    dn7 = kd1 * vitKa * Xa * Vol_liver
    dn10 = dn7
    dXa = n12 / Vol_liver - kd1*vitKa*Xa

    #box E
    dthr = (-ek1*thr*fbng) + (ek_1*ES) + (ek2*ES) - (k_thrombin_deg*thr) + (n10/Vol_blood)
    dfbng = (-ek1*thr*fbng) + (ek_1*ES) + (n8/Vol_blood)
    dES = (ek1*thr*fbng) - (ek_1*ES) - (ek2*ES)
    dn14 = (ek2*ES*Vol_blood)

    #box B
    dVitKi = (n3/Vol_liver) - (k_vitKa_exc*vitKi) - (n5/Vol_liver) + (n7/Vol_liver)

    #box C
    dn5 = (ck*vitKi*Vol_liver)
    dn6 = (ck*vitKi*Vol_liver)

    #box F
    #dTF = (-5/10**6)*R
    dTF = -tf1*0.1
    dX = (n11/Vol_blood) - (R/Vol_blood)
    dn12 = R
    dfbr = (n14/Vol_blood) - k_fibrin_deg*fbr

    return np.array([dn5-n5, dn6-n6, dn7-n7, dn10-n10, dn12-n12, dn14-n14,
                     dVitKa, dVitKi, dX, dXa, dthr, 
                     dfbng, dES, dTF, dfbr])

output2 = odeint(odefunc2, y0, tspan)

vitKaResult2 = output2[:,6]
XResult2 = output2[:,8]
XaResult2 = output2[:,9]
thrResult2 = output2[:,10]
fbrResult2 = output2[:,14]
n6Result2 = output2[:,1]

print(tspan)
print(fbrResult)
print(fbrResult2)

print("max:")
print(np.max(fbrResult2))

fig, ax1 = plt.subplots()
#ax1.plot(tspan, n6, label="n6")
ax1.plot(tspan, fbrResult, label="Healthy", color='steelblue')
ax1.plot(tspan, fbrResult2, label="Diseased", color='crimson')
plt.axhline(y = 18313.43, color = 'k', linestyle = 'dashed', label="Desired Level of Fibrin for Clotting")
#plt.axhline(y = 9000, color = 'k', linestyle = '-', label="Maintenance Level")
#ax1.plot(tspan, cES, label="cES")
#ax1.plot(tspan, R2, label="R2")
#ax1.plot(tspan, L, label="L")
plt.title("Decreased Fibrin in Diseased State Yields Increased Clotting Time")
#plt.ylim(0,20000)
#plt.xlim
plt.xlabel("Time (sec)")
plt.ylabel("Fibrin Concentration (nM)")
plt.legend(loc = "best")
plt.show()

thrombin_max = np.max(thrResult)
thrombin_maxD = np.max(thrResult2)
print("thrombin max healthy: {:3e}".format(thrombin_max))
print("thrombin max diseased: {:3e}".format(thrombin_maxD))

fig, ax1 = plt.subplots()
#ax1.plot(tspan, n6, label="n6")
#ax1.plot(tspan, XaResult, label="Xa Healthy")
#ax1.plot(tspan, XaResult2, label="Xa Diseased")
ax1.plot(tspan, thrResult, label="Healthy", color='crimson')
ax1.plot(tspan, thrResult2, label="Diseased")
#plt.axhline(y = 16200, color = 'k', linestyle = 'dashed', label="Desired Level of Fibrin for Clotting")
#plt.axhline(y = 9000, color = 'k', linestyle = '-', label="Maintenance Level")
#ax1.plot(tspan, cES, label="cES")
#ax1.plot(tspan, R2, label="R2")
#ax1.plot(tspan, L, label="L")
plt.title("Vitamin K Deficiency Leads to Decreased Levels of Thrombin")
#plt.ylim(0,20000)
plt.xlim(0, 60)
plt.xlabel("Time (sec)")
plt.ylabel("Thrombin Concentration (nM)")
plt.legend(loc = "best")
plt.show()

"""
fig, ax1 = plt.subplots()
#ax1.plot(tspan, n6, label="n6")
ax1.plot(tspan, n6Result, label="Healthy")
ax1.plot(tspan, n6Result2, label="Diseased")
plt.axhline(y = 1.541*10, color = 'k', linestyle = 'dashed', label="Steady State")
#plt.axhline(y = 16200, color = 'k', linestyle = 'dashed', label="Desired Level of Fibrin for Clotting")
#plt.axhline(y = 9000, color = 'k', linestyle = '-', label="Maintenance Level")
#ax1.plot(tspan, cES, label="cES")
#ax1.plot(tspan, R2, label="R2")
#ax1.plot(tspan, L, label="L")
plt.title("Vitamin K (Active) in Stream 6 Decreases in Disease State")
#plt.ylim(0,20000)
plt.xlim(0, 20)
plt.xlabel("Time (sec)")
plt.ylabel("Flow rate in stream 6 (nmol/sec)")
plt.legend(loc = "best")
plt.show()"""



X = ['0.4','0.7','1.0','1.3']
healthy = [8.78, 14.23, 20.22, 30.2]
diseased = [13.63, 23.28, 34.59, 50.8]
  
X_axis = np.arange(len(X))
  
plt.bar(X_axis - 0.2, healthy, 0.4, label = 'Healthy', color='steelblue')
plt.bar(X_axis + 0.2, diseased, 0.4, label = 'Diseased', color='darkslategray')
  
plt.xticks(X_axis, X)
plt.xlabel("Length of Wound (cm)")
plt.ylabel("Clotting Time (sec)")
plt.title("Clotting Time Increases as Injury Size Increases (Width = 0.5cm")
plt.legend()
plt.show()

vitKSpaced = np.linspace(0, 1.45, 100)
thr = 0.68955*vitKSpaced
fbr = 6206.67*vitKSpaced

fig, ax = plt.subplots(2)

#ax1.plot(tspan, n6, label="n6")
ax[0].plot(vitKSpaced, thr, label="Thrombin", color="crimson")
ax[1].plot(vitKSpaced, fbr, label="Fibrin", color="darkslategray")
ax[0].invert_xaxis()
ax[1].invert_xaxis()
ax[0].set(ylabel="Thrombin Concentration (nmol)", 
          title='Fibrin and Thrombin Concentrations Increase as Vitamin K Active Increases')
ax[1].set(xlabel='Vitamin K Active Concentration (nmol)', ylabel='Fibrin Concentration (nmol)')
plt.show()

differences = [4.82, 9.10, 14.4, 20.6]
plt.bar(X, differences, color="darkslategray")
plt.title('Disparity Between Clotting Time in Healthy\nand Diseased Increases as Injury Size Increases')
plt.xlabel('Injury Width (cm)')
plt.ylabel('Difference in Clotting Time (%)')
plt.show()