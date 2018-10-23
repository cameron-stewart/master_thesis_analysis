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

plotwidth = 4.5 

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

########## Stress Plots ##################################################################################
fig, axes = plt.subplots(4,2,figsize=(6.,9.5))
plt.subplots_adjust(hspace=0.3)
box_props = dict(boxstyle='round', facecolor=colors[2], alpha=0.8)
for i in range(len(plots)):
    plot = plots[i]
    y_sim = y[(L/2-npoints):-(L/2-npoints)]/float(L)
    y_ana = yfit/float(L)

    #rescale xx component
    stressxx_ana = Axx(yfit, alpha[plot], lambdap[plot], etap, popt[plot,0])
    stressxx_sim = stress[plot,(L/2-npoints):-(L/2-npoints),0]
    stressxx_ana = stressxx_ana-np.amin(stressxx_sim)
    stressxx_sim = stressxx_sim-np.amin(stressxx_sim)
    stressxx_ana = 0.8*stressxx_ana/np.amax(stressxx_sim)+0.1
    stressxx_sim = 0.8*stressxx_sim/np.amax(stressxx_sim)+0.1
     
    # rescale yy component
    stressyy_ana = Ayy(yfit, alpha[plot], lambdap[plot], etap, popt[plot,1])
    stressyy_sim = stress[plot,(L/2-npoints):-(L/2-npoints),2]
    stressyy_ana = stressyy_ana-np.amax(stressyy_sim)
    stressyy_sim = stressyy_sim-np.amax(stressyy_sim)
    stressyy_ana = -0.8*stressyy_ana/np.amin(stressyy_sim)-0.1
    stressyy_sim = -0.8*stressyy_sim/np.amin(stressyy_sim)-0.1
    
  
    #plot xx
    axes[i,0].plot(y_sim, stressxx_sim, 'o', color=colors[0], label='Simulated')
    axes[i,0].plot(y_ana, stressxx_ana, color=colors[3], label='Analytical')
    axes[i,0].set_xticks(np.arange(-.04, 0.06, 0.04))
    axes[i,0].set_yticks(np.arange(0.0, 1.1, 0.5))
    axes[i,0].set_ylim(0.0,1.0)
    if i!=3:
        axes[i,0].set_xticklabels([])
    
    axes[i,0].spines["top"].set_visible(False)    
    axes[i,0].spines["right"].set_visible(False)    
    #axes[i,0].spines["bottom"].set_visible(False)    
    #axes[i,0].spines["left"].set_visible(False)    
    axes[i,0].get_xaxis().tick_bottom()
    axes[i,0].get_yaxis().tick_left()
    axes[i,0].text(0.04,0.85,r'$\epsilon =\ $' + '%.2f' %round(epsilon[plot],2),transform=axes[i,0].transAxes, fontsize=15, bbox=box_props)

    
    #plot yy
    axes[i,1].plot(y_sim, stressyy_sim, 'o', color=colors[1], label='Simulated')
    axes[i,1].plot(y_ana, stressyy_ana, color=colors[3], label='Analytical')
    axes[i,1].set_xticks(np.arange(-.04, 0.06, 0.04))
    axes[i,1].set_yticks(np.arange(-1.0, 0.1, 0.5))
    axes[i,1].set_ylim(-1.0,0.0)
    if i!=3:
        axes[i,1].set_xticklabels([])
    
    axes[i,1].spines["top"].set_visible(False)    
    axes[i,1].spines["left"].set_visible(False)
    #axes[i,1].spines["bottom"].set_visible(False)    
    #axes[i,1].spines["right"].set_visible(False)    
    axes[i,1].get_xaxis().tick_bottom()
    axes[i,1].get_yaxis().tick_right()

fig.text(0.515,0.04, r'$y/L$',ha='center')
fig.text(0.515,0.92,'Incompressible',ha='center',fontsize=22)
fig.text(0.02,0.5,r'$\tau_{xx}(0,y)\ (\mathrm{AU})$', va='center', rotation='vertical')
fig.text(1.0,0.5,r'$\tau_{yy}(0,y)\ (\mathrm{AU})$', va='center', rotation=-90)

plt.savefig('/home/cstewart/thesis/plots/roll_stress_inc.pgf', bbox_inches='tight')
