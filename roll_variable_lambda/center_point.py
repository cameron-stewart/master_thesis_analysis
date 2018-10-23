import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import spline
import sys
from matplotlib.legend_handler import HandlerBase

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

########## Center Point Stress ###########################################################################

# lambda arrays
x1_ana = epsilon[:17]
x1_sim = epsilon
y1_ana = 2*etap*alpha[:17]/(1-2*epsilon[:17])
y1_sim = stress[:,L/2+1,0]
x2_ana = epsilon
x2_sim = epsilon
y2_ana = -2*etap*alpha/(1+2*epsilon)
y2_sim = stress[:,L/2+1,2]

# force arrays
x1_ana_f = np.load('x1_ana_f.npy')
x1_sim_f = np.load('x1_sim_f.npy')
y1_ana_f = np.load('y1_ana_f.npy')
y1_sim_f = np.load('y1_sim_f.npy')
x2_ana_f = np.load('x2_ana_f.npy')
x2_sim_f = np.load('x2_sim_f.npy')
y2_ana_f = np.load('y2_ana_f.npy')
y2_sim_f = np.load('y2_sim_f.npy')

# interpolation
x1_new = np.linspace(x1_ana.min(), x1_ana.max(),300)
y1_new = spline(x1_ana,y1_ana,x1_new)
x1_new_f = np.linspace(x1_ana_f.min(), x1_ana_f.max(),300)
y1_new_f = spline(x1_ana_f,y1_ana_f,x1_new_f)

# normalize
norm = y1_sim_f.max()*1.05

y1_new = y1_new/norm
y1_new_f = y1_new_f/norm
y1_sim = y1_sim/norm

y2_ana = y2_ana/norm
y2_ana_f = y2_ana_f/norm
y2_sim = y2_sim/norm
y2_sim_f = y2_sim_f/norm

y1_sim_f = y1_sim_f/norm

fig,axes = plt.subplots(2,1,figsize=(6,8),sharex=True)
ax1 = axes[0]
ax2 = axes[1]

ax1.plot(x1_new, y1_new,color=colors[0], label='Analytical', linewidth=2)
ax1.plot(x1_sim, y1_sim, 'o', color=colors[0], label='Simulated')
ax1.plot(x1_new_f, y1_new_f,color=colors[1], label='Analytical', linewidth=2)
ax1.plot(x1_sim_f, y1_sim_f, 'o', color=colors[1], label='Simulated')
ax1.axvline(0.25,color='black',linewidth=1.0, linestyle='--')
ax1.axvline(0.5,color='black',linewidth=1.0, linestyle='--')
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()
ax1.set_ybound(0.0,1.0)
ax1.set_ylabel(r'$\tau_{xx}(0,0)/\tau_{xx}^\mathrm{max}$')
ax1.text(0.14,-0.07,'Smooth')
ax1.text(0.35,-0.07,'Cusp')
ax1.text(0.55,-0.07,'Singularity')

ax2.plot(x2_ana, y2_ana,color=colors[0], label='Analytical', linewidth=2)
ax2.plot(x2_sim, y2_sim, '^', color=colors[0], label='Simulated')
ax2.plot(x2_ana_f, y2_ana_f,color=colors[1], label='Analytical', linewidth=2)
ax2.plot(x2_sim_f, y2_sim_f, '^', color=colors[1], label='Simulated')
ax2.set_xlabel(r'$\epsilon=\alpha\lambda_p$')
ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.get_xaxis().tick_bottom()
ax2.get_yaxis().tick_left()
ax2.set_ylabel(r'$\tau_{yy}(0,0)/\tau_{xx}^\mathrm{max}$')

class AnyObjectHandler(HandlerBase):
   def create_artists(self, legend, orig_handle,
                      x0, y0, width, height, fontsize, trans):
       l1 = plt.Line2D([x0+8.,x0+8.], [0.7*height,0.7*height],
                       linestyle=orig_handle[1], color=orig_handle[0], marker='o')
       l2 = plt.Line2D([x0,y0+width], [0.0*height,0.0*height],
                       linestyle=orig_handle[1], color=orig_handle[0], linewidth=2)
       l3 = plt.Line2D([x0+20.,x0+20.], [0.7*height,0.7*height],
                       linestyle=orig_handle[1], color=orig_handle[0], marker='^')
       return [l1, l2, l3]

ax1.legend([(colors[1],'-'), (colors[0],"-")], ['Varying' + '$\ f$', 'Varying' + '$\ \lambda_p$'],
           handler_map={tuple: AnyObjectHandler()}, loc=2)

plt.savefig('/home/cstewart/thesis/plots/roll_center.pgf', bbox_inches='tight')

