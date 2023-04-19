import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp

# %% Initial Values
#reaction for vitmain k reduction (box c)
#rate constants

#degradation rate constants
k_vitKa_exc = (1.868*10**-3)*3600

#Reaction Rate
R = (1.232*10**2) #1/day

#Streams
n1=  266   #nmol/day
n3 = (((1.693*10**-3))/(2))*3600 #nmol/sec
n5 = (1.541*10)     #nmol/sec
n7 = (1.541*10)     #nmol/sec

#known values for steady state
vitKa = 1.45                #nmol/L
vitKi = vitKa                #nmol/L
Vol_blood = 5               #Liters
Vol_liver = Vol_blood*0.125 #Liters


# %% Calculations
tspan = np.linspace(0, 20, 10000) #timpoint for each second 
y0 = np.array([vitKi])

def odefunc(y0, t): 
    vitKi = y0[0]

    #box B
    dVitKi = (n3/Vol_liver) - (k_vitKa_exc*vitKi) - (n5/Vol_liver) + (n7/Vol_liver)

    return np.array([dVitKi])

output = odeint(odefunc, y0, tspan)

vitKiResult = output[:,0]

fig, ax1 = plt.subplots()
#ax1.plot(tspan, n6, label="n6")
ax1.plot(tspan, vitKiResult, label="vitKi")
#ax1.plot(tspan, cES, label="cES")
#ax1.plot(tspan, R2, label="R2")
#ax1.plot(tspan, L, label="L")
#plt.title("Concentrations of Double Receptor-Ligand\nReaction vs. Time")
#plt.xlim(0,120)
plt.xlabel("Time (sec)")
plt.ylabel("Concentration (M)")
plt.legend(loc = "best")
plt.show()
