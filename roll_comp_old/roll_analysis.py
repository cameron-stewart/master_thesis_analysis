import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys

def Axx(y, alpha, lambdap, etap, C):
    epsilon = alpha*lambdap
    return C*np.power(np.absolute(y),(1-2*epsilon)/(epsilon)) + 2*etap*alpha/(1 - 2*epsilon)

def Ayy(y, alpha, lambdap, etap, C):
    epsilon = alpha*lambdap
    return C*np.power(np.absolute(y),2.0+1/epsilon) - 2*etap*alpha/(1+2*epsilon)

lambdap = np.array([8000,10000,12000,13000,14000,15000,16000,17000,18000,19000,20000,21000,22000,23000,24000,26000,28000,32000])
lambdap = np.array([10000,18000,24000,28000])
L = 213
etap = 0.075
npoints = 8
nruns = lambdap.size
stress = np.zeros((nruns,L,3))
vel = np.zeros((nruns,L))
alpha = np.zeros(nruns)
epsilon = np.zeros(nruns)
exponent = np.zeros(nruns)
popt = np.zeros((nruns,2))
pcov = np.zeros((nruns,2))
tyy_pi = np.zeros(nruns)
y = np.arange(-(L/2),L/2+1)

for j in range(nruns):
    with open('stress' + str(lambdap[j]) + '.csv','rb') as f:
        reader = csv.reader(f)
        i = 0
        for row in reader:
            stress[j,i,0] = float(np.array(row)[0])
            stress[j,i,1] = float(np.array(row)[1])
            stress[j,i,2] = float(np.array(row)[4])
            i += 1

    with open('vel' + str(lambdap[j]) + '.csv','rb') as f:
        reader = csv.reader(f)
        i = 0
        for row in reader:
            vel[j,i] = float(np.array(row)[1])
            i += 1

    alpha[j] = (vel[j,L/2-1]-vel[j,L/2+1])/2.0
    epsilon[j] = alpha[j]*lambdap[j]
    exponent[j] = ((1-2*epsilon[j])/epsilon[j])

    popt[j,0],pcov[j,0] = curve_fit(lambda y, C: Axx(y, alpha[j], lambdap[j], etap, C), 
            y[L/2-npoints:L/2], stress[j,L/2-npoints:L/2,0])

    popt[j,1],pcov[j,1] = curve_fit(lambda y, C: Ayy(y, alpha[j], lambdap[j], etap, C), 
            y[L/2-npoints:L/2], stress[j,L/2-npoints:L/2,2])

yfit = np.linspace(-npoints, npoints, 100)
wi = lambdap*1e-3/L

########## Stress Plots ##################################################################################
for plot in range(nruns):  
    plt.figure(figsize = (8/1.5,6/1.5))
    plt.plot(y[(L/2-npoints):-(L/2-npoints)]/float(L), stress[plot,(L/2-npoints):-(L/2-npoints),0], 'o', label='Simulated')
    plt.plot(yfit/float(L), Axx(yfit, alpha[plot], lambdap[plot], etap, popt[plot,0]), 'k', label='Analytical')
    plt.xlabel('$y/L$', size=20)
    plt.ylabel(r'$\tau_{xx}(0,y)\ (\mathrm{LU})$', size=20)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend(numpoints=1, loc=0)
    plt.xticks(np.arange(-.04,0.06, 0.02))
    plt.title(r'$\epsilon =$' + '%.2f' %round(epsilon[plot],2))
    plt.savefig('/home/cstewart/thesis/figures/comp_stressxx_' + str(lambdap[plot]) + '_roll.eps', 
            format='eps', dpi=1200, bbox_inches='tight')

for plot in range(nruns):  
    plt.figure(figsize = (8/1.5,6/1.5))
    plt.plot(y[(L/2-npoints):-(L/2-npoints)]/float(L), stress[plot,(L/2-npoints):-(L/2-npoints),2], 'o', label='Simulated')
    plt.plot(yfit/float(L), Ayy(yfit, alpha[plot], lambdap[plot], etap, popt[plot,1]), 'k', label='Analytical')
    plt.xlabel('$y/L$', size=20)
    plt.ylabel(r'$\tau_{yy}(0,y)\ (\mathrm{LU})$', size=20)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xticks(np.arange(-.04,0.06, 0.02))
    plt.legend(numpoints=1, loc=0)
    plt.title(r'$\epsilon =$' + '%.2f' %round(epsilon[plot],2))
    plt.savefig('/home/cstewart/thesis/figures/comp_stressyy_' + str(lambdap[plot]) + '_roll.eps', 
            format='eps', dpi=1200, bbox_inches='tight')
