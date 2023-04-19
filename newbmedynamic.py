import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp

# %% Initial Values
#reaction for vitmain k reduction (box c)
#rate constants
ck =
#reaction for thrombin activation (box d)
#rate constants
dk =
#reaction for injury site (box f)
#rate constants
fk =
#reaction for fibrin gen (box e)
#rate constants
ek =
#degradation rate constants
k_thrombin_deg =     # 1/sec
k_fibrin_deg =      # 1/sec
xa_deg =          #1/sec

#Streams
n1=  0.0030787037   #nmol/sec
n2 =    #nmol/sec
n3 =    #nmol/sec
n4 =    #nmol/sec
n5 =      #nmol/sec
n6 =      #nmol/sec
n7 =      #nmol/sec
n8 =     #nmol/sec
n9 =     #nmol/sec
n10 =    #nmol/sec
n11 =    #nmol/sec
n12 =     #nmol/sec
n13 =     #nmol/sec
n14 =    #nmol/sec
n15 =    #nmol/sec

#known values for steady state
vitKa = 1.45                #nmol/L
Xa = 170                    #nmol/L
proth = 1400                #nmol/L
fbng = 9000*10**3           #nmol/L
thr = 1*10**3               #nmol/L
tf1 = 1.6*10**-5            #nmol/L
X = 160                     #nmol/L
Xi = 160                    #nmol/L
vitKi = 1.45                #nmol/L
fbr = fbng                  #nmol/L
Vol_blood = 5               #Liters
Vol_liver = Vol_blood*0.125 #Liters
cES = 0
ddES = cES
eES = cES
VKORC1 = 9.857*10**-3       #nmol/L     #MAYBE CHANGE !!!!!!!!!
Xprime = X

#rate constants
kc2 = 3*10**6 #1/sec
kd22 = 57 #1/sec
kf1 = 0.025 #L/sec*nnmol
ke2= 84 #1/sec