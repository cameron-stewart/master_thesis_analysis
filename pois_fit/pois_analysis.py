
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
force = 1.48e-5
l = 1801.8
L = 30
eta_s=0.5
eta_p=0.5
umax = force*L*L/(8.*(eta_s+eta_p))


def parabola(Ly, F, etas, etap, y):
    return -(F/2)*(y*y - Ly*y)/(etas+etap)

def vgrad(Ly, F, etas, etap, y):
    return -(F/2)*(2*y-Ly)/(etas+etap)

def stressxx_ana(Ly, F, etas, etap, lambdap, y):
    return 2*lambdap*etap*(vgrad(Ly, F, etas, etap, y))**2

def stressxy_ana(Ly, F, etas, etap, lambdap, y):
    return etap*vgrad(Ly, F, etas, etap, y)


y_fit = np.linspace(0,L,100)
vel = np.zeros(L)
stressxx = np.zeros(L)
stressxy = np.zeros(L)
stressyy = np.zeros(L)
y = np.zeros(L)

with open('vel.csv','rb') as f:
    reader = csv.reader(f)
    i = 0
    for row in reader:
        vel[i] = float(np.array(row)[0])
        y[i] = float(np.array(row)[6]) - 0.5
        i += 1

with open('stress.csv','rb') as f:
    reader = csv.reader(f)
    i = 0
    for row in reader:
        stressxx[i] = float(np.array(row)[0])
        stressxy[i] = float(np.array(row)[1])
        stressyy[i] = float(np.array(row)[4])
        i += 1


fig,axes = plt.subplots(2,2,figsize=(plotwidth,0.8*plotwidth))
ax1 = axes[0,0]
ax2 = axes[0,1]
ax3 = axes[1,0]
ax4 = axes[1,1]

ax1.plot(y/L, vel/umax, 'o', label='Simulated', marker='o', markersize=5, color=colors[0])
ax1.plot(y_fit/L,parabola(L,force,eta_s,eta_p,y_fit)/umax, color='black', label='Analytical', linewidth=1.2)
#ax1.spines["top"].set_visible(False)
#ax1.spines["right"].set_visible(False)
#ax1.get_xaxis().tick_bottom()
#ax1.get_yaxis().tick_left()
ax1.set_ylabel(r'$u_x/u_x^\mathrm{max}$')
ax1.set_xticks([0,.25,.5,.75,1])
ax1.set_yticks([0.0,0.3,0.6,0.9,1.2])
ax1.tick_params(labelbottom=False)   

ax2.plot(y/L, stressxx/force, 'o', label='Simulated', marker='o', markersize=5, color=colors[1])
ax2.plot(y_fit/L,stressxx_ana(L,force,eta_s,eta_p,l,y_fit)/force, color='black', label='Analytical', linewidth=1.2)
#ax2.spines["top"].set_visible(False)
#ax2.spines["left"].set_visible(False)
#ax2.get_xaxis().tick_bottom()
ax2.get_yaxis().tick_right()
ax2.set_ylabel(r'$\tau_{xx}/(f_x/\delta x^2)$')
ax2.yaxis.set_label_position("right")
ax2.set_xticks([0,.25,.5,.75,1])
ax2.set_yticks([0,1.5,3,4.5,6])
ax2.tick_params(labelbottom=False,labelleft=False,labelright=True)   

ax3.plot(y/L, stressxy/force, 'o', label='Simulated', marker='o', markersize=5, color=colors[2])
ax3.plot(y_fit/L,stressxy_ana(L,force,eta_s,eta_p,l,y_fit)/force, color='black', label='Analytical', linewidth=1.2)
#ax3.spines["top"].set_visible(False)
#ax3.spines["right"].set_visible(False)
#ax3.get_xaxis().tick_bottom()
#ax3.get_yaxis().tick_left()
ax3.set_ylabel(r'$\tau_{xy}/(f_x/\delta x^2)$')
ax3.set_xticks([0,.25,.5,.75,1])
ax3.set_yticks([-8,-4,0,4,8])

ax4.plot(y/L, stressyy/force, 'o', label='Simulated', marker='o', markersize=5, color=colors[4])
ax4.plot(y_fit/L,np.zeros(100)/force, color='black', label='Analytical', linewidth=1.2)
#ax4.spines["top"].set_visible(False)
#ax4.spines["left"].set_visible(False)
#ax4.get_xaxis().tick_bottom()
ax4.get_yaxis().tick_right()
ax4.set_ylim(-1,1)
ax4.set_xticks([0,.25,.5,.75,1])
ax4.set_ylabel(r'$\tau_{yy}/(f_x/\delta x^2)$')
ax4.yaxis.set_label_position("right")
ax4.tick_params(labelleft=False,labelright=True)   


fig.text(0.515,0.01, r'$y/L$',ha='center')

plt.savefig('/home/cstewart/thesis/plots/pois_fit.pgf', bbox_inches = 'tight')

