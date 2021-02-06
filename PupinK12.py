#Semi-empirical equations from "TR-11 AERODYNAMIC DRAG OF MODEL ROCKET by Gerald M. Gregorek#
#Using standard atmosphere

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Main paramaters
L_bt=2.7 #[m]-Body tube length
L_n=0.3 #[m]-Nose length
D=0.120 #[m]- Main diameter

#Stabilizator paramaters
n=0.5#Tamper ratio [/] 0<n<1
n_s=4 #[/]-Number of stabilizators(fins)
l_s=0.15 #[m]-Root cord stabilizator
h_s=0.15 #[m]-Stabilizator span
t_c=0.3 #[/]-Airfoil thicknes ratio (max thic/root chod)

#Air constants
kappa=1.4
R=287

######################End of parameter imput################################

##Standar Atmospehere##
h=np.linspace(0,4200,14) #Height [m]
rho=1.225*(1-(h/44300))**(4.25588) # Air density[kg/m^3] 
T=288.15-0.0065*h #Temprature [K]
ni=1.458*10**(-6)*T**(1.5)/(T+110.4) #Dynamic Viscosity [Pa*s]
M=np.linspace(0.05,0.7,14) #Mach Number [/]
V=M*np.sqrt(kappa*R*T)# Speed [m/s]
#Rocket Geometry

L=L_bt+L_n #Whole lennght of the rocket [m]
L_D=L/D #Rocket Aspect Ratio [/]

S_ref=D**2*np.pi/4 #Reference cross section [m^2]
S_bt=D*np.pi*L_bt #Body tube area [m^2]
S_f=2*1.02*(l_s*(1+n))/2*h_s*n_s #Total Fin area [m^2]
S_wet=D*np.pi*((L_n**2+D**2)/2+L_bt+D/4)+(l_s*(1+n))/2*h_s*n#Wetted Area [m^2]


Re_bt=rho*V*L/ni #Body tube Reynolds Number
Re_f=rho*V*l_s/ni #Fin Reynolds Number

Cf_bt=0.455/(np.log10(Re_bt))**2.58 #Skin drag coefficient
Cf_f=0.455/(np.log10(Re_f))**2.58 #Skin drag coefficient

Cdn_Cdbt=1.02*Cf_bt*(1+(1.5/((L_D)**(3/2))))*S_wet/S_ref #Body tube and nose drag coeff
Cdb=0.029/(np.sqrt(Cdn_Cdbt)) #Base drag coefficient
Cdf=2*Cf_f*(1+2*t_c)*S_f/S_ref #Fin drag coefficient
Cdint=2*Cf_f*l_s/S_ref*D/2*n_s #Interirence drag coeffitient

Cd0theoretical=Cdn_Cdbt+Cdb+Cdf+Cdint #Total zero drag coeffitient
Cd0=Cd0theoretical*1.3 #Zero drag coeff
F=0.5*rho*V**2*S_ref*Cd0/10 #Drag force [daN]


#Rocket stability (No boattail or difference in body tube diameter
#Nose
Cnalfan=2 #Normal force coefficient of nose (Gereeral is 2)
Xn=0.466*L_n #Distance between the nose tip and cp of the nose (See reference [2])

#Fin
l=np.sqrt((1/2*l_s*(1-n))**2+h_s**2) #Mid cord line
Cnalfaf=2*n_s*(h_s/D)**2/(1+np.sqrt(1+(2*l/(l_s*(1+n)))**2)) #Normal force coeeficient of fin
Kfb=1+D/2/(h_s+D/2) #Interference factor (for n= 3 or 4)
Cnalfaf=Cnalfaf*Kfb #Corrected normal force coefficient of fin
m=l_s*(1-n)
xf=L-l_s #Distance between tip of nose and start of root cord of fins
x_f=xf+m/2*(1*n)/(1+n)+1/6*(l_s*(1+n-l_s*n/(1+n))) #Fin centere of presure

Cnalfa=Cnalfan+Cnalfaf #Normal force coefficinet of the rocket
Xcp=(Cnalfan*Xn+Cnalfaf*x_f)/Cnalfa # Center of presure mesured from the nose tip [m]
Xcprel=Xcp/L #Relative center of presure

Xcgrel=np.linspace(0,1,21) #Relative center of mass
Xcg=Xcgrel*L #Center of mass [m]
Stab=(Xcprel-Xcgrel)*L/D #Stability [cal]
Xcg_Stab=np.dstack((Xcg,Stab)) #Creates a matrix

#Graphs
plt.figure(1)
plt.plot(M,Cd0,linewidth=3)
plt.grid()
plt.title('$Drag\ Coefficient$',fontsize=15)
plt.xlabel('$Mach\ Number$',fontsize=15)
plt.ylabel('$C_{D_{0}}$',fontsize=15)


plt.figure(2)
plt.plot(M,F,linewidth=3)
plt.grid()
plt.title('$Drag\ Force$',fontsize=15)
plt.xlabel('$Mach\ Number$',fontsize=15)
plt.ylabel('$F_D[daN]$',fontsize=15)

plt.figure(3)
plt.plot(Stab,Xcg,linewidth=3)
plt.grid()
plt.title('$Change\ of\ stability$',fontsize=15)
plt.xlabel('$Stability\ [cal]$',fontsize=15)
plt.ylabel('$CG\ position [m]$',fontsize=15)
plt.show()


#Saving resoults of the  drag calculation into cvs file#
resultsdrag=pd.DataFrame({"M":M, "Cdn_Cdbt":Cdn_Cdbt,"Cdf":Cdf,"Cdint":Cdint, "Cd0":Cd0,"F":F})
resultsdrag.to_csv("resultsdrag.cvs",index=False)
#Savig results of the stability calculation into csv file#
resultsstability=pd.DataFrame({"Xcg [m]":Xcg, "Stab [cal]":Stab, "Fin span [m]":h_s, "Root cord[m]":l_s,"Tamper ratio [/]":n, "Number of fins [/]":n_s})
resultsstability.to_csv("resultsstability.csv",index=False)


#Saving all results
results=pd.concat([resultsdrag,resultsstability])
results.to_csv("results.csv",index=False)