import math
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.signal import argrelextrema
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

plotwidth = 6 
leng = 10000
betas = np.array([1,3,5,7,9])
betas = np.array([9,7,5,3,1])
nruns = betas.size

t1 = np.zeros(nruns)
vmax = np.zeros(nruns)
t2 = np.zeros(nruns)
v2 = np.zeros(nruns)
vel = np.load('numpy.save.npy')
x = np.arange(leng)

vel = vel/vel[-1,-1]
vel_norm = vel - 1.0

for i in range(nruns):
    vmax[i] = np.amax(vel[:,i])
    t1[i] = np.argmax(vel[:,i])
    for t in range(len(vel[:,i])):
        if ((vel_norm[t,i] < vel_norm[int(t1[i]),i]/math.e) & ( t > int(t1[i]))):
            t2[i] = t- t1[i]
            v2[i] = vel[t,i]
            break 

plt.figure(figsize = (plotwidth,plotwidth*0.75))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

for i in range(nruns):
    plt.plot(x,vel[:,i],color=colors[i],label=(r'$\beta =\ $'+"0." + str(betas[i])), linewidth=2.0)
    plt.plot(t1,vmax,'d',color='k', markersize=4.0)
    plt.plot(t2+t1,v2,'d',color=colors[5], markersize=4.0)

"""
plt.plot(x,vel[:,-1],color=colors[nruns-1],label=(r'$\beta =\ $'+"0." + str(betas[-1])), linewidth=2.0)
plt.plot(t1[-1],vmax[-1],'d',color='k', markersize=4.0)
plt.plot(t2[-1]+t1[-1],v2[-1],'d',color=colors[5], markersize=4.0)
plt.plot(x,vel[:,-2],color=colors[nruns-2],label=(r'$\beta =\ $'+"0." + str(betas[-2])), linewidth=2.0)
plt.plot(t1[-2],vmax[-2],'d',color='k', markersize=4.0)
plt.plot(t2[-2]+t1[-2],v2[-2],'d',color=colors[5], markersize=4.0)
"""
plt.ylim(0.8,1.4)
plt.xlim(0,8000)
plt.legend(frameon=False)
plt.xlabel("Time Steps" + r'$/\delta t$')
plt.ylabel('Relative Velocity' + r'\ $u_x/u_x^\mathrm{ss}$')
"""
ax2=plt.axes([.40,.6,.2,.2])
ax2.plot(betas/10.0,t1,'d-',color='k', markersize=4.0)
ax2.set_ylabel(r'$t_1/\delta t$', color='k')
#ax2.set_yticks([300,500,700])
ax2.tick_params('y',colors='k')
ax2.set_xlabel(r'$\beta$')

a=ax2.twinx()
a.plot(betas/10.0,t2, '-d', color=colors[5], markersize=4.0)
a.set_yticks([800,1600,2400])
a.tick_params('y', colors=colors[5])
a.set_ylabel(r'$t_2/\delta t$', color=colors[5])
"""
#plt.savefig("/home/cstewart/thesis/plots/overshoot_beta.pgf", bbox_inches="tight")
plt.savefig("overshoot_beta.pdf", bbox_inches='tight')
