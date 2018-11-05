import math
import csv
import numpy as np
import matplotlib.pyplot as plt
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

plotwidth = 6 
leng = 10000
lambdas = np.array([1000,3000,5000,7000,9000])
nruns = lambdas.size

t1 = np.zeros(nruns)
vmax = np.zeros(nruns)
t2 = np.zeros(nruns)
v2 = np.zeros(nruns)
#vel = np.zeros((leng,nruns))
vel = np.load("numpy.save.npy")
x = np.arange(0,leng)

vel = vel/vel[-1,0]
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
    plt.plot(vel[:,i],color=colors[i],label=(r'$\lambda_p =\ $' + str(lambdas[i])),linewidth=2.0)
    plt.plot(t1,vmax,'d',color='k', markersize=4.0)
    plt.plot(t2+t1,v2,'d',color=colors[5], markersize=4.0)

plt.ylim(0.7,1.4)
plt.legend(frameon=False)
plt.xlabel("Time Steps" + r'$\ (\delta t)$')
plt.ylabel('Relative Velocity' + r'$\ u_x/u_x^\mathrm{ss}$')

ax2=plt.axes([.42,.2,.2,.2])
ax2.plot(lambdas,t1,'d-',color='k', markersize=4.0)
ax2.set_ylabel(r'$t_1 (\delta t)$', color='k')
ax2.set_yticks([300,500,700])
ax2.tick_params('y',colors='k')
ax2.set_xlabel(r'$\lambda_p\ (\delta t)$')

a=ax2.twinx()
a.plot(lambdas,t2, '-d', color=colors[5], markersize=4.0)
a.set_yticks([1500,4500,7500])
a.tick_params('y', colors=colors[5])
a.set_ylabel(r'$t_2 (\delta t)$', color=colors[5])


plt.savefig("/home/cstewart/thesis/plots/overshoot_lambda.pgf", bbox_inches="tight")
#plt.savefig("overshoot_lambda.pdf", bbox_inches="tight")






"""
plt.plot(betas/10.,vmax, 'o')
plt.axis([0,1,1,1.8])
plt.plot(betas/10.,tmax, 'o')
plt.axis([0,1,300,500])
"""
