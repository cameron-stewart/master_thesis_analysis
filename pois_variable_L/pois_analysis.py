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

plotwidth = 4.5 

def parabola(Ly, F, etas, etap, y):
    return -(F/2)*(y*y - Ly*y)/(etas+etap)

def vgrad(Ly, F, etas, etap, y):
    return -(F/2)*(2*y-Ly)/(etas+etap)

def stressxx(Ly, F, etas, etap, lambdap, y):
    return 2*lambdap*etap*(vgrad(Ly, F, etas, etap, y))**2

def stressxy(Ly, F, etas, etap, lambdap, y):
    return etap*vgrad(Ly, F, etas, etap, y)


y30 = np.linspace(0,30,100)
y300 = np.linspace(0,300,100)
vel1 = np.zeros(30)
vel2 = np.zeros(300)
stressxx1 = np.zeros(30)
stressxy1 = np.zeros(30)
stressyy1 = np.zeros(30)
stressxx2 = np.zeros(300)
stressxy2 = np.zeros(300)
stressxx3 = np.zeros(30)
stressyy3 = np.zeros(30)
y1 = np.zeros(30)
y2 = np.zeros(300)

with open('velocity_30.csv','rb') as f:
    reader = csv.reader(f)
    i = 0
    for row in reader:
        vel1[i] = float(np.array(row)[0])
        y1[i] = float(np.array(row)[6]) - 0.5
        i += 1

with open('velocity_300.csv','rb') as f:
    reader = csv.reader(f)
    i = 0
    for row in reader:
        vel2[i] = float(np.array(row)[1])
        y2[i] = float(np.array(row)[7]) - 0.5
        i += 1

with open('stress_30.csv','rb') as f:
    reader = csv.reader(f)
    i = 0
    for row in reader:
        stressxx1[i] = float(np.array(row)[0])
        stressxy1[i] = float(np.array(row)[1])
        stressyy1[i] = float(np.array(row)[4])
        i += 1

with open('stress_300.csv','rb') as f:
    reader = csv.reader(f)
    i = 0
    for row in reader:
        stressxx2[i] = float(np.array(row)[0])
        stressxy2[i] = float(np.array(row)[1])
        i += 1
with open('stress_01.csv','rb') as f:
    reader = csv.reader(f)
    i = 0
    for row in reader:
        stressxx3[i] = float(np.array(row)[0])
        stressyy3[i] = float(np.array(row)[4])
        i += 1


fig1 = plt.figure(figsize = (plotwidth,plotwidth*.75))
ax1 = plt.subplot(111)
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()
plt.plot(y1, vel1, 'o', label='Simulated', color=colors[5])
plt.plot(y30,parabola(30,1e-3,0.5,0.5,y30), 'k', label='Analytical')
plt.xlabel(r'$y\ (\mathrm{LU})$', size = 20)
plt.ylabel(r'$u_x\ (\mathrm{LU})$', size = 20)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(numpoints = 1, loc = 0)
plt.savefig('/home/cstewart/thesis/plots/pois_velocity_30.pgf', bbox_inches = 'tight')

fig2 = plt.figure(figsize = (plotwidth,plotwidth*.75))
ax2 = plt.subplot(111)
ax2.get_xaxis().tick_bottom()
ax2.get_yaxis().tick_left()
plt.plot(y2[4::10], vel2[4::10], 'o', label='Simulated', color=colors[2])
plt.plot(y300,parabola(300,1e-6,0.5,0.5,y300), 'k', label='Analytical')
plt.xlabel(r'$y\ (\mathrm{LU})$', size = 20)
plt.ylabel(r'$u_x\ (\mathrm{LU})$', size = 20)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(numpoints = 1, loc = 0)
plt.savefig('/home/cstewart/thesis/plots/pois_velocity_300.pgf', bbox_inches = 'tight')

fig3 = plt.figure(figsize = (plotwidth,plotwidth*.75))
ax3 = plt.subplot(111)
ax3.get_xaxis().tick_bottom()
ax3.get_yaxis().tick_left()
plt.plot(y1, stressxy1, 'o', label='Simulated', color=colors[5])
plt.plot(y30,stressxy(30,1e-3,0.5,0.5, 300, y30), 'k', label='Analytical')
plt.xlabel(r'$y\ (\mathrm{LU})$', size = 20)
plt.ylabel(r'$\tau_{xy}\ (\mathrm{LU})$', size = 20)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(numpoints = 1, loc = 0)
plt.savefig('/home/cstewart/thesis/plots/pois_stressxy_30.pgf', bbox_inches = 'tight')

fig4 = plt.figure(figsize = (plotwidth,plotwidth*.75))
ax4 = plt.subplot(111)
ax4.get_xaxis().tick_bottom()
ax4.get_yaxis().tick_left()
plt.plot(y2[4::10], stressxy2[4::10], 'o', label='Simulated', color=colors[2])
plt.plot(y300,stressxy(300,1e-6,0.5,0.5,30000,y300), 'k', label='Analytical')
plt.xlabel(r'$y\ (\mathrm{LU})$', size = 20)
plt.ylabel(r'$\tau_{xy}\ (\mathrm{LU})$', size = 20)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(numpoints = 1, loc = 0)
plt.savefig('/home/cstewart/thesis/plots/pois_stressxy_300.pgf', bbox_inches = 'tight')

fig5 = plt.figure(figsize = (plotwidth,plotwidth*.75))
ax5 = plt.subplot(111)
ax5.get_xaxis().tick_bottom()
ax5.get_yaxis().tick_left()
plt.plot(y1, stressxx1, 'o', label='Simulated, Wi $= 1.125$', color=colors[0])
plt.plot(y30,stressxx(30,1e-3,0.5,0.5, 300, y30), 'k')
plt.plot(y1, stressxx3, 'o', label='Simulated, Wi $= 0.1125$', color=colors[1])
plt.plot(y30, stressxx(30,1e-3,0.5,0.5,30,y30), 'k', label = 'Analytical')
plt.xlabel(r'$y\ (\mathrm{LU})$', size = 20)
plt.ylabel(r'$\tau_{xx}\ (\mathrm{LU})$', size = 20)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(numpoints = 1, loc = 0)
plt.savefig('/home/cstewart/thesis/plots/pois_stressxx_wi.pgf', bbox_inches = 'tight')

fig6 = plt.figure(figsize = (plotwidth,plotwidth*.75))
ax6 = plt.subplot(111)
ax6.get_xaxis().tick_bottom()
ax6.get_yaxis().tick_left()
plt.plot(y1, stressyy1, 'o', label='Simulated, Wi $= 1.125$', color=colors[0])
plt.plot(y30,np.zeros(100), 'k')
plt.plot(y1, stressyy3, 'o', label='Simulated, Wi $= 0.1125$', color=colors[1])
plt.plot(y30, np.zeros(100), 'k', label = 'Analytical')
plt.xlabel(r'$y\ (\mathrm{LU})$', size = 20)
plt.ylabel(r'$\tau_{yy}\ (\mathrm{LU})$', size = 20)
plt.axis([0,30,-5e-15,5e-15])
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(numpoints = 1, loc = 0)
plt.savefig('/home/cstewart/thesis/plots/pois_stressyy_wi.pgf', bbox_inches = 'tight')


