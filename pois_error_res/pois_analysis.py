
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

def parabola(Ly, F, etas, etap, y):
    return -(F/2)*(y*y - Ly*y)/(etas+etap)

def vgrad(Ly, F, etas, etap, y):
    return -(F/2)*(2*y-Ly)/(etas+etap)

def stressxx_ana(Ly, F, etas, etap, lambdap, y):
    return 2*lambdap*etap*(vgrad(Ly, F, etas, etap, y))**2

def stressxy_ana(Ly, F, etas, etap, lambdap, y):
    return etap*vgrad(Ly, F, etas, etap, y)


nruns = 4
force = np.full(nruns,1.48e-5)
lambdap = np.full(nruns,1801.8)
L = np.full(nruns,30,dtype=int)
umax = np.empty(nruns)
eta_s=0.5
eta_p=0.5

for i in range(nruns):
    L[i] = L[i]*(i+1)
    force[i] = force[i]/(i+1)**3
    lambdap[i] = lambdap[i]*(i+1)**2
    umax[i] = force[i]*L[i]**2/(8.*(eta_s+eta_p))

y_fit = np.zeros((100,nruns))
vel = np.zeros((L[-1],nruns))
stressxx = np.zeros((L[-1],nruns))
stressxy = np.zeros((L[-1],nruns))
stressyy = np.zeros((L[-1],nruns))
y = np.zeros((L[-1],nruns))

for j in range(nruns):
    l = L[j]
    y_fit[:,j] = np.linspace(0,l,100)
    with open('vel' + str(l) + '.csv','rb') as f:
        reader = csv.reader(f)
        i = 0
        for row in reader:
            vel[i,j] = float(np.array(row)[0])
            y[i,j] = float(np.array(row)[6]) - 0.5
            i += 1

    with open('stress' + str(l) + '.csv','rb') as f:
        reader = csv.reader(f)
        i = 0
        for row in reader:
            stressxx[i,j] = float(np.array(row)[0])
            stressxy[i,j] = float(np.array(row)[1])
            stressyy[i,j] = float(np.array(row)[4])
            i += 1


fig,axes = plt.subplots(2,2,figsize=(plotwidth,0.8*plotwidth))
ax1 = axes[0,0]
ax2 = axes[0,1]
ax3 = axes[1,0]
ax4 = axes[1,1]

run = 0

ax1.plot(y[:L[run],run]/L[run], vel[:L[run],run]/umax[run], 'o', label='Simulated', marker='o', markersize=5, color=colors[0])
ax1.plot(y_fit[:,run]/L[run],parabola(L[run],force[run],eta_s,eta_p,y_fit[:,run])/umax[run], color='black', label='Analytical', linewidth=1.2)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()
ax1.set_ylabel(r'$u_x/u_x^\mathrm{max}$')
ax1.set_xticks([0,.25,.5,.75,1])
ax1.set_yticks([0.0,0.3,0.6,0.9,1.2])

ax2.plot(y[:L[run],run]/L[run], stressxx[:L[run],run]/(force[run]), 'o', label='Simulated', marker='o', markersize=5, color=colors[1])
ax2.plot(y_fit[:,run]/L[run],stressxx_ana(L[run],force[run],eta_s,eta_p,lambdap[run],y_fit[:,run])/(force[run]), color='black', label='Analytical', linewidth=1.2)
ax2.spines["top"].set_visible(False)
ax2.spines["left"].set_visible(False)
ax2.get_xaxis().tick_bottom()
ax2.get_yaxis().tick_right()
ax2.set_ylabel(r'$\tau_{xx}\ (f_x/\delta x^2)$')
ax2.yaxis.set_label_position("right")
ax2.set_xticks([0,.25,.5,.75,1])
#ax2.set_yticks([0,2,4,6,8])

ax3.plot(y[:L[run],run]/L[run], stressxy[:L[run],run]/force[run], 'o', label='Simulated', marker='o', markersize=5, color=colors[2])
ax3.plot(y_fit[:,run]/L[run],stressxy_ana(L[run],force[run],eta_s,eta_p,lambdap[run],y_fit[:,run])/force[run], color='black', label='Analytical', linewidth=1.2)
ax3.spines["top"].set_visible(False)
ax3.spines["right"].set_visible(False)
ax3.get_xaxis().tick_bottom()
ax3.get_yaxis().tick_left()
ax3.set_ylabel(r'$\tau_{xy}\ (f_x/\delta x^2)$')
ax3.set_xticks([0,.25,.5,.75,1])
#ax3.set_yticks([-8,-4,0,4,8])

ax4.plot(y[:L[run],run]/L[run], stressyy[:L[run],run]/force[run], 'o', label='Simulated', marker='o', markersize=5, color=colors[2])
ax4.plot(y_fit[:,run]/L[run],np.zeros(100), color='black', label='Analytical', linewidth=1.2)
ax4.spines["top"].set_visible(False)
ax4.spines["left"].set_visible(False)
ax4.get_xaxis().tick_bottom()
ax4.get_yaxis().tick_right()
ax4.set_ylim(-1,1)
ax4.set_xticks([0,.25,.5,.75,1])
ax4.set_ylabel(r'$\tau_{yy}\ (f_x/\delta x^2)$')
ax4.yaxis.set_label_position("right")

fig.text(0.515,0.01, r'$y/L$',ha='center')

plt.savefig('/home/cstewart/thesis/plots/pois_test.pdf', bbox_inches = 'tight')

