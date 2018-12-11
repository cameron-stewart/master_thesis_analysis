import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys

import matplotlib as mpl
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['text.usetex'] = True
mpl.rcParams['pgf.rcfonts'] = False
mpl.rcParams['pgf.preamble'] = [r'\usepackage{bm}']
mpl.rcParams['pgf.texsystem'] = 'pdflatex'
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{bm}',r'\usepackage[T1]{fontenc}']


colors = [(0,107,164),(255,128,14),(171,171,171),(89,89,89),(95,158,209),(200,82,0),(137,137,137),(162,200,236),(255,188,121),(207,207,207)]

for i in range(len(colors)):
    r,g,b = colors[i]
    colors[i] = (r/255.,g/255.,b/255.)

plotwidth = 5

def Axx(y, alpha, lambdap, etap, C):
    epsilon = alpha*lambdap
    return C*np.power(np.absolute(y),(1-2*epsilon)/(epsilon)) + 2*etap*alpha/(1 - 2*epsilon)

def Ayy(y, alpha, lambdap, etap, C):
    epsilon = alpha*lambdap
    return C*np.power(np.absolute(y),2.0+1/epsilon) - 2*etap*alpha/(1+2*epsilon)

lambdap = np.array([8000,10000,12000,13000,14000,15000,16000,17000,18000,19000,20000,21000,22000,23000,24000,26000,28000,32000,36000,
    40000,44000])
plots = np.array([1,8,14,17])
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

########## Parameters ####################################################################################


fig, ax1 = plt.subplots(figsize = (plotwidth,plotwidth*0.75))
alpha = alpha*1e5
ax1.plot(lambdap,alpha,'o-', label='Simulated',color=colors[0])
ax1.set_xlabel(r'$\lambda_p/\delta t$')
ax1.set_ylabel(r'$\alpha\times 10^5/(1/\delta t)$', color=colors[0])
ax1.set_xbound(5000,45000)
ax1.set_xticks(np.arange(10000,45001,10000))
ax1.tick_params('y',colors=colors[0])

ax2=ax1.twinx()
ax2.plot(lambdap,epsilon,'o-', label='Simulated',color=colors[1])
ax2.set_xbound(5000,45000)
ax2.set_xticks(np.arange(10000,45001,10000))
ax2.set_ylabel(r'$\epsilon = \alpha\lambda_p$', color=colors[1])
ax2.tick_params('y',colors=colors[1])

ax3=ax1.twiny()
ax3.set_xlim(ax1.get_xlim()[0]/(1000.*float(L)),ax1.get_xlim()[1]/(1000.*float(L)))
ax3.set_xlabel(r'$\mathrm{Wi}$')
plt.savefig('/home/cstewart/thesis/plots/roll_parameters.pgf',bbox_inches='tight')

