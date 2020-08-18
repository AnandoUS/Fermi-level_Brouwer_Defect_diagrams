import matplotlib.pyplot as plt
import numpy as np
from fdint import * 
from scipy.integrate import quad

########################CONSTANTS####################################################
kB = 1.38064852 * 10**-23		#Boltzman Constant SI unit
e = 1.6*10**-19				#Charge on electron SI unit
h = 6.626*10**-34			#Planks Constant SI unit
pi = np.pi				#pi
m_e = 9.10938356 * 10**-31		#mass of an electron in SI unit

T=773

mp = 0.409	         		#Seebeck mass of Valence band
mn = 0.542				#Seebeck mass of Conduction band 

v = 132.44/10**24			#Volume of a primitive unit-cell in cm-3
Eg = 0.428				#Band-gap in eV

Efs = 0.1				#Starting Fermi-level
Efe = 0.3				#Ending Fermi-level
npts = 2000				#Number of points between starting and ending Ef 

c = (kB*T)/e				parameter to convert fermi-level to reduced ferm-level

########### Initiating arrays for ferm-level and defect energies #####################

Efermi = np.array([x for x in np.linspace(Efs,Efe,npts)])
ChemPot = np.array([x for x in np.linspace(0,0.823,npts)])

V_Mg = np.array([x for x in np.linspace(1.616,0.789,npts)])
Mg_i = np.array([x for x in np.linspace(0.039,0.862,npts)])
Sb_Mg = np.array([x for x in np.linspace(2.569,0.512,npts)])

p = np.array([((4*pi)*(((2*mp*m_e*kB*T)/(h**2))**1.5)*(fdk(0.5,(-Ef)/(c))))/(10**6) for Ef in np.linspace(Efs,Efe,npts)])		#hole carrier concentration in cm-3
n = np.array([((4*pi)*(((2*mn*m_e*kB*T)/(h**2))**1.5)*(fdk(0.5,(Ef-Eg)/(c))))/(10**6) for Ef in np.linspace(Efs,Efe,npts)])		#electron carrier concentration in cm-3

############################ Initiating empty arrays #################################
EF=[]
VMg = []
Mgi = []
SbMg = []
ne = []
ph = []

for i in range(npts):
	DefV_Mg = np.array([(3/v)*np.exp(-(V_Mg[i]-2*Ef)/(c)) for Ef in np.linspace(Efs,Efe,npts)])
	DefMg_i = np.array([(1/v)*np.exp(-(Mg_i[i]+2*Ef)/(c)) for Ef in np.linspace(Efs,Efe,npts)])
	DefSb_Mg = np.array([(1/v)*np.exp(-(Sb_Mg[i]+1*Ef)/(c)) for Ef in np.linspace(Efs,Efe,npts)])

	k = np.argmin(np.abs(p-n + 2*DefMg_i - 2*DefV_Mg + 1*DefSb_Mg))						#Numerically soving for charge neutrality
	
	EF.append(Efermi[k])

	VMg.append((3/v)*np.exp(-(V_Mg[i]-2*Efermi[k])/(c)))
	Mgi.append((1/v)*np.exp(-(Mg_i[i]+2*Efermi[k])/(c)))
	SbMg.append((1/v)*np.exp(-(Sb_Mg[i]+1*Efermi[k])/(c)))
	ne.append(((4*pi)*(((2*mn*m_e*kB*T)/(h**2))**1.5)*(fdk(0.5,(Efermi[k]-Eg)/(c))))/(10**6))
	ph.append(((4*pi)*(((2*mp*m_e*kB*T)/(h**2))**1.5)*(fdk(0.5,(-Efermi[k])/(c))))/(10**6))

plt.plot(ChemPot,VMg)
plt.plot(ChemPot,Mgi)
plt.plot(ChemPot,SbMg)
plt.plot(ChemPot,ne,linewidth=4.5)
plt.plot(ChemPot,ph,linewidth=4.5)

plt.yscale('log')
plt.ylim(10**15,10**19)

plt.xlabel('Mg Atomic Chemical Potential (eV)', fontsize = 16)
plt.ylabel('Defect Concentrations (cm$^{-3}$)', fontsize = 16)

plt.show()

plt.plot(ChemPot,EF)
plt.show()