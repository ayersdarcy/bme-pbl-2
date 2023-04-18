import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp

# %% Initial Values
#reaction for vitmain k reduction (box c)
#rate constants
cK1 = 2 #1/nM * s
cK_1 = 3 #1/s
cK2 = 5 #1/s
#reaction for thrombin activation (box d)
#rate constants
dK1 = 8  #1/nM*s
ddK1 = 5     #1/nM*s
ddK_1 = 4         #1/s
ddK2 = 6           #1/s
#reaction for injury site (box f)
#rate constants
fK1 = 3   #1/nM*s
#reaction for fibrin gen (box e)
#rate constants
eK1 = 4   #1/nM*s
eK_1 = 6      #1/s
eK2= 2         #1/s
#degradation rate constants
k_thrombin_deg = 2    # 1/sec
k_fibrin_deg = 3     # 1/sec
xa_deg = 6           #1/sec

#Streams
n1= 3    #nmol/sec
n2 = 3   #nmol/sec
n3 = 3   #nmol/sec
n4 = 3   #nmol/sec
n5 = 4     #nmol/sec
n6 = 4     #nmol/sec
n7 = 4     #nmol/sec
n8 = 3    #nmol/sec
n9 = 2    #nmol/sec
n10 = 3   #nmol/sec
n11 = 6   #nmol/sec
n12 = 4    #nmol/sec
n13 = 4    #nmol/sec
n14 = 2   #nmol/sec
n15 = 4    #nmol/sec
n16 = 2   #nmol/sec

#known values for steady state
vitKa = 1.45                #nmol/L
Xa = 12                    #nmol/L
proth = 14                #nmol/L
fbng = 3           #nmol/L
thr = 3               #nmol/L
tf1 = 5            #nmol/L
X = 3                     #nmol/L
Xi = 5                    #nmol/L
vitKi = 1.45                #nmol/L
fbr = fbng                  #nmol/L
Vol_blood = 5               #Liters
Vol_liver = Vol_blood*0.125 #Liters
cES = 0
ddES = cES
eES = cES
VKORC1 = 3       #nmol/L     #MAYBE CHANGE !!!!!!!!!
Xprime = X

#rate constants
kc2 = 6 #1/sec
kd22 = 57 #1/sec
kf1 = 0.025 #L/sec*nnmol
ke2= 84 #1/sec


# %% Calculations 
tspan = np.linspace(0, 10, 100) #timpoint for each minute 
y0 = np.array([n4, n5, n6, n7, n11, n13, n16, vitKa, vitKi, VKORC1, 
      cES, Xprime, Xa, proth, ddES, tf1, X, fbr, thr, fbng, eES])

print("started")

def odefunc(y0, t):

    #value definition
    n1 = 3
    n2 = 3
    n3 = 3
    n4 = y0[0]
    n5 = y0[1]
    n6 = y0[2]
    n7 = y0[3]
    n8 = 3
    n9 = 2
    #n10 not used for calculations 
    n11 = y0[4]
    n12 = 4
    n13 = y0[5]
    #n14 and n15 aren't used for these calculation
    n16 = y0[6]

    vitKa = y0[7]
    vitKi = y0[8]
    VKORC1 = y0[9]
    cES = y0[10]
    Xprime = y0[11]
    Xa = y0[12]
    proth = y0[13]
    ddES = y0[14]
    tf1 = y0[15]
    X = y0[16]
    fbr = y0[17]
    thr = y0[18]
    fbng = y0[19]
    eES = y0[20]

    #box a & b (Liver Production)
    dn4 = 0.55 * (n3 + n7)
    dn5 = 0.45 * (n3 + n7)
    #n8, n12, and n9 don't change because they are all generation terms

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

    #box e (capillaries)
    deThr = -eK1*thr*fbng + eK_1*eES + eK2*eES - k_thrombin_deg*thr + n11/Vol_blood
    deFbng = -eK1*thr*fbng + eK_1*eES + n9/Vol_blood
    deES = eK1*thr*fbng - eK1*eES - eK2*eES
        #deFbr = eK2*eES - n16/Vol_blood
    dn16 = (eK2*eES)*Vol_blood
    
    return np.array([dn4, dn5, dn6, dn7, dn11, dn13, dn16, 
            dVitKa, dVitKi, dVKORC1, dcES, dXprime,
            dXa + ddXa, ddproth, ddES, dTF, dfX, dffbr,
            deThr, deFbng, deES])


output = odeint(odefunc, y0, tspan)

#fbr = output[:,17]
#vitKi = output[:,8]
#n6 = output[:,2]
#cES = output[:,10]
fbr = output[:,17]
#CConcentration = output["y"][1]

fig, ax1 = plt.subplots()
#ax1.plot(tspan, n6, label="n6")
ax1.plot(tspan, fbr, label="vitKi")
#ax1.plot(tspan, cES, label="cES")
#ax1.plot(tspan, R2, label="R2")
#ax1.plot(tspan, L, label="L")
#plt.title("Concentrations of Double Receptor-Ligand\nReaction vs. Time")
plt.xlabel("Time (sec)")
plt.ylabel("Concentration (M)")
plt.legend(loc = "best")
plt.show()

print("Done")