import matplotlib.pyplot as plt
import numpy as np
from fdint import * 
from scipy.integrate import quad

##################### PbTe data from James's paper###################################
data=np.loadtxt('Carrier_vs_Pb_Vacancy_Conc.txt')
data1=np.loadtxt('Nominal_Carrier_vs_Pb_Vacancy_Conc.txt')

########################CONSTANTS####################################################
kB = 1.38064852 * 10**-23		#Boltzman Constant SI unit
e = 1.6*10**-19				#Charge on electron SI unit
h = 6.626*10**-34			#Planks Constant SI unit
pi = np.pi				#pi
m_e = 9.10938356 * 10**-31		#mass of an electron in SI unit

T=973

mp = 0.468	
mn = 0.26*(T/300)**0.5				#Seebeck mass of Conduction band 	

v = 67/10**24				#Volume of a primitive unit-cell in cm-3
Eg = 0.4				#Band-gap in eV
Efs = -0.5				#Starting Fermi-level
Efe = 1					#Ending Fermi-level
npts = 1500				#Number of points between starting and ending Ef 

Efermi = np.array([x for x in np.linspace(Efs,Efe,npts)])

##################################################################################### High Temperature  #################################################################################
a = (((kB*T)/e)/Eg) 
c = ((kB*T)/e)

def integrand(x, l,m,n, Ef):
	return ((((1+2*a*x)**2+2)**(l/2))*((x+a*x**2)**m)*((x)**n))*(np.exp(x-(Ef/c))/(1+np.exp(x-(Ef/c)))**2)

F__2_1_0 = np.array([quad(integrand, 0, 350, args=(-2,1,0,ef-Eg))[0] for ef in np.linspace(Efs,Efe,npts)])
F__4_05_0 = np.array([quad(integrand, 0, 350, args=(-4,0.5,0,Ef-Eg))[0] for Ef in np.linspace(Efs,Efe,npts)])

###################################   Carrier Concentrations     #######################################################
p = np.array([((16/3)*pi*(((2*mp*m_e*kB*T)/(h**2))**1.5)*(fdk(0,(-Ef)/(c))**2/fdk(-0.5,(-Ef)/(c))))/(10**6) for Ef in np.linspace(Efs,Efe,npts)])	#electron carrier concentration in cm-3
n = np.array((((4)*pi*(((2*mn*m_e*kB*T)/(h**2))**1.5))*((F__2_1_0)**2/(F__4_05_0)))/(10**6))

###################################  Calculating Defect Energies and fermi-level   ################################################
a=0.82				
Def = np.array([(a-2*Ef) for Ef in np.linspace(Efs,Efe,npts)])												#Defect Energy at VBM
conc = (1/v)*np.exp(-(Def)/(c))																#defect concentration in cm-3

cconc=data[:,0]
carry=(p-n)/10**18
exper1=[]

for i in cconc:
	ja = np.argmin(np.abs(i-carry))
	exper1.append(Efermi[ja])

exper1 = np.array(exper1)

c=1.487*10**4
VPb = np.array(data[:,1]*c)
defc=np.array([(1/v)*np.exp(-(en)/((kB*T)/e))/10**18 for en in np.linspace(0.3,0.8,npts)])
EN = np.array([en for en in np.linspace(0.3,0.8,npts)])
exper2 = [] 

for i in VPb:
	ja = np.argmin(np.abs(i-defc))
	exper2.append(EN[ja])

exper2 = np.array(exper2)

############################### Defect energy Diagram ################################
a1=0.82

dE1 = [a1-2*x for x in np.linspace(-1,0.5,5)]
x1 = [x for x in np.linspace(-1,0.5,5)]

plt.plot(x1, dE1,c = 'steelblue', linewidth=1.5)
plt.plot([0,0],[0,1.5],c='steelblue',linewidth=1.5)
plt.plot([0.4,0.4],[0,1.5],c='steelblue',linewidth=1.5)

plt.scatter(exper1,exper2,s=80)

plt.xlim(-0.1,0.5)
plt.ylim(0,1)

plt.rcParams['figure.figsize'] = 8, 8

plt.xticks(fontsize= 14)
plt.yticks(fontsize= 14)

plt.xlabel('Fermi-level (eV)', fontsize= 16)
plt.ylabel('Defect energy (eV)', fontsize =15)

plt.show()