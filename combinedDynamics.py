import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp

# %% Initial Values
#reaction for vitmain k reduction (box c)
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
tspan = np.linspace(0, 60, 10000) #timpoint for each second 
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

vitKiResult = output[:,14]


# CHANGE TO DISEASED DYNAMIC 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp

# %% Initial Values
#reaction for vitmain k reduction (box c)
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
tspan = np.linspace(0, 60, 10000) #timpoint for each second 
y0 = np.array([n5, n6, n7, n10, n12, n14,
                vitKa, vitKi, X, Xa, thr, 
                fbng, ES, tf1, fbr])

def odefunc2(y0, t): 
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

vitKiResult2 = output2[:,14]

fig, ax1 = plt.subplots()
#ax1.plot(tspan, n6, label="n6")
ax1.plot(tspan, vitKiResult, label="vitKi")
ax1.plot(tspan, vitKiResult2, label="diseased")
plt.axhline(y = 16200, color = 'r', linestyle = 'dashed')
#ax1.plot(tspan, cES, label="cES")
#ax1.plot(tspan, R2, label="R2")
#ax1.plot(tspan, L, label="L")
#plt.title("Concentrations of Double Receptor-Ligand\nReaction vs. Time")
#plt.ylim(0,20000)
#plt.xlim
plt.xlabel("Time (sec)")
plt.ylabel("Concentration (M)")
plt.legend(loc = "best")
plt.show()