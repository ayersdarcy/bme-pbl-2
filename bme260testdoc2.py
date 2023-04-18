import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp

# %% Initial Values
#reaction for vitmain k reduction (box c)
#rate constants
cK1 = 1000 #1/nM * s
cK_1 = 10 #1/s
cK2 = 5 #1/s
#reaction for thrombin activation (box d)
#rate constants
dK1 = 2.207*10**-6  #1/nM*s
ddK1 = 1*10**-1     #1/nM*s
ddK_1 = 103         #1/s
ddK2 = 57           #1/s
#reaction for injury site (box f)
#rate constants
fK1 = 2.5*10**5 / (10**9)   #1/nM*s
#reaction for fibrin gen (box e)
#rate constants
eK1 = 1*10**2   #1/nM*s
eK_1 = 636      #1/s
eK2= 84         #1/s
#degradation rate constants
k_thrombin_deg = 5.652*10**2    # 1/sec
k_fibrin_deg = 5.185*10**-3     # 1/sec
xa_deg = 3.012*10**-6           #1/sec

#Streams
n1= 3.079*10**-3    #nmol/sec
n2 = 1.385*10**-3   #nmol/sec
n3 = 1.693*10**-3   #nmol/sec
n4 = 1.693*10**-3   #nmol/sec
n5 = 3.2*10**-4     #nmol/sec
n6 = 3.2*10**-4     #nmol/sec
n7 = 3.2*10**-4     #nmol/sec
n8 = 2.826*10**3    #nmol/sec
n9 = 2.333*10**2    #nmol/sec
n10 = 2.826*10**3   #nmol/sec
n11 = 2.826*10**6   #nmol/sec
n12 = .32       #nmol/sec
n13 = 3.2*10**-4    #nmol/sec
n14 = 2.333*10**2   #nmol/sec
n15 = 3.2*10**-4    #nmol/sec
n16 = 2.333*10**2   #nmol/sec

#known values for steady state
vitKa = 1.45                #nmol/L
Xa = 170                    #nmol/L
proth = 1400                #nmol/L
fbng = 9000                 #nmol/L
thr = 1*10**3               #nmol/L
tf1 = 16                   #nmol/L
X = 160                     #nmol/L
Xi = 160                    #nmol/L
vitKi = 1.45                #nmol/L
fbr = fbng                  #nmol/L
Vol_blood = 5               #Liters
Vol_liver = Vol_blood*0.125 #Liters
cES = 0
ddES = cES
eES = cES
VKORC1 = 5*10**-3       #nmol/L     #MAYBE CHANGE !!!!!!!!!
Xprime = X

#rate constants
kc2 = 3*10**6 #1/sec
kd22 = 57 #1/sec
kf1 = 4.4*10**-20 #L/sec*nnmol
ke2= 84 #1/sec


# %% Calculations 
tspan = np.linspace(0, 100, 1000) #timpoint for each minute 
y0 = np.array([X, Xprime, tf1, fbr, vitKi, VKORC1, cES, n6, Xprime, Xa, n7, vitKa, proth, ddES, n11 ])

print("started")

def odefunc(y0, t):

    X = y0[0]
    n13 = y0[1]
    tf1 = y0[2]
    fbr = y0[3]
    vitKi = y0[4]
    VKORC1 = y0[5]
    cES = y0[6]
    n6 = y0[7]
    Xprime = y0[8]
    Xa = y0[9]
    n7 = y0[10]
    vitKa =y0[11]
    proth = y0[12]
    ddES = y0[13]
    n11 = y0[14]

    #box c (vit k reduction)
    dVitKi = -cK1*vitKi*VKORC1 + cK_1*cES + n5/Vol_liver
        #n5 = Vol_liver * (cK1*vitKi*VKORC1 - cK_1*cES) #not solving for this anymore, solved in box B
    dVKORC1 = -cK1*vitKi*VKORC1 + cK_1*cES + cK2*cES 
    dcES = cK1*vitKi*VKORC1 - cK_1*cES - cK2*cES
        #dcVitKa_activation = cK2*cES - n6/Vol_liver
    dn6 = Vol_liver * cK2*cES

    #box d (thrombin act)
    dXprime = -dK1* Xprime * vitKa + n13/Vol_liver
    dXa = dK1*Xprime*vitKa
        #dVitKi_X = dK1*Xprime*vitKi - n7/Vol_liver
    dn7 = Vol_liver*(dK1*Xprime*vitKa)
    dVitKa = -dK1*Xprime*vitKi + n6/Vol_liver
    
    ddXa = -ddK1*Xa*proth + ddK_1*ddES + ddK2*ddES - k_thrombin_deg*Xa
    ddproth = -ddK1*Xa*proth + ddK_1*ddES + n8/Vol_liver 
    #n8 won't change bc it's proth generation
    ddES = ddK1*Xa*proth - ddK_1*ddES -ddK2*ddES
        #ddthr = ddK2*ddES - n11/Vol_liver
    dn11 = (ddK2*ddES)*Vol_liver

    #box f (injury site)
    dTF = -fK1*tf1*X 
    dfX = -fK1*tf1*X + n12/Vol_blood
    dffbr = n16/Vol_blood - k_fibrin_deg*fbr #added this fibrin equation?
    #n12 is x_gen so it is not changing
        #dfXprime = fK1*tf1*X - n13/Vol_blood
    dn13 = (fK1*tf1*X)*Vol_blood

    return np.array([dfX, dn13, dTF, dffbr, 
                     dVitKi, dVKORC1, dcES, dn6,
                     dXprime, dn7, dVitKa, ddXa+dXa, ddproth, ddES, dn11])

output = odeint(odefunc, y0, tspan)

dXprime = output[:,8]

fig, ax1 = plt.subplots()
ax1.plot(tspan, dXprime, label="dffbr")
plt.title("dprime vs Time")
plt.xlabel("Time (sec)")
#plt.ylim([0,300])
plt.ylabel("Concentration (M)")
plt.legend(loc= "best")
plt.show()

# %%
