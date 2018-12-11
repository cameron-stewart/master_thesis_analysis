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

plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelleft=False) # labels along the bottom edge are off

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
for i in range(nruns):
    plt.plot(x,vel[:,i],color=colors[i],label=(r'$\beta =\ $'+"0." + str(betas[i])), linewidth=2.0)



plt.savefig("/home/cstewart/thesis/inkscape/title.pdf", bbox_inches="tight")
#plt.savefig("overshoot_beta.pdf", bbox_inches='tight')
