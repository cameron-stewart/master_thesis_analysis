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

force = np.array([5,7,8,10,12,13,15,17,18,20,22,23])
lambdap = 18000
L = 213
etap = 0.075
npoints = 8
nruns = force.size
velmax = np.zeros(nruns)
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
    with open('stress' + str(force[j]) + '.csv','rb') as f:
        reader = csv.reader(f)
        i = 0
        for row in reader:
            stress[j,i,0] = float(np.array(row)[0])
            stress[j,i,1] = float(np.array(row)[1])
            stress[j,i,2] = float(np.array(row)[4])
            i += 1

    with open('vel' + str(force[j]) + '.csv','rb') as f:
        reader = csv.reader(f)
        i = 0
        for row in reader:
            vel[j,i] = float(np.array(row)[1])
            i += 1
        velmax[j] = vel[j,:].max()

    alpha[j] = (vel[j,L/2-1]-vel[j,L/2+1])/2.0
    epsilon[j] = alpha[j]*lambdap
    exponent[j] = ((1-2*epsilon[j])/epsilon[j])

    popt[j,0],pcov[j,0] = curve_fit(lambda y, C: Axx(y, alpha[j], lambdap, etap, C), 
            y[L/2-npoints:L/2], stress[j,L/2-npoints:L/2,0])

    popt[j,1],pcov[j,1] = curve_fit(lambda y, C: Ayy(y, alpha[j], lambdap, etap, C), 
            y[L/2-npoints:L/2], stress[j,L/2-npoints:L/2,2])

yfit = np.linspace(-npoints, npoints, 100)
print(velmax*float(L)/0.15)

### ######### Center Point Stress ###########################################################################

x1_ana = epsilon[:7]
x1_sim = epsilon
y1_ana = 2*etap*alpha[:7]/(1-2*epsilon[:7])
y1_sim = stress[:,L/2+1,0]

x2_ana = epsilon
x2_sim = epsilon
y2_ana = -2*etap*alpha/(1+2*epsilon)
y2_sim = stress[:,L/2+1,2]

np.save('x1_ana_f',x1_ana)
np.save('x1_sim_f',x1_sim)
np.save('y1_ana_f',y1_ana)
np.save('y1_sim_f',y1_sim)
np.save('x2_ana_f',x2_ana)
np.save('x2_sim_f',x2_sim)
np.save('y2_ana_f',y2_ana)
np.save('y2_sim_f',y2_sim)

"""
plt.figure(figsize = (8/1.5,6/1.5))
plt.plot(epsilon,-2*etap*alpha/(1+2*epsilon),'k', label='Analytical')
plt.plot(epsilon, stress[:,L/2+1,2], 'o', label='Simulated')
plt.xlabel(r'$\epsilon=\alpha\lambda_p$', size = 20)
plt.ylabel(r'$\tau_{yy}(0,0)$', size = 20)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(numpoints=1, loc = 0)
plt.savefig('/home/cstewart/thesis/plots/roll_stressyy_center_alpha.eps', format='eps', dpi=1200, bbox_inches='tight')

plt.figure(figsize = (8/1.5,6/1.5))
plt.plot(epsilon,np.zeros(nruns),'k', label='Analytical')
plt.plot(epsilon, stress[:,L/2+1,1], 'o', label='Simulated')
plt.xlabel(r'$\epsilon=\alpha\lambda_p$', size = 20)
plt.ylabel(r'$\tau_{xy}(0,0)$', size = 20)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(numpoints=1, loc = 0)
plt.axis([0.1,0.7,-1e-15,1e-15])
plt.savefig('/home/cstewart/thesis/plots/roll_stressxy_center_alpha.eps', format='eps', dpi=1200, bbox_inches='tight')

plt.figure(figsize = (8/1.5,6/1.5))
plt.plot(epsilon[:7],2*etap*alpha[:7]/(1-2*epsilon[:7]),'k', label='Analytical')
plt.plot(epsilon, stress[:,L/2+1,0], 'o', label='Simulated')
plt.xlabel(r'$\epsilon=\alpha\lambda_p$', size = 20)
plt.ylabel(r'$\tau_{xx}(0,0)$', size = 20)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(numpoints=1, loc = 0)
plt.savefig('/home/cstewart/thesis/plots/roll_stressxx_center_alpha.eps', format='eps', dpi=1200, bbox_inches='tight')
"""
