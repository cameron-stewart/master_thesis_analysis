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


def parabola(Ly, F, etas, etap, y):
    return -(F/2)*(y*y - Ly*y)/(etas+etap)

def vgrad(Ly, F, etas, etap, y):
    return -(F/2)*(2*y-Ly)/(etas+etap)

def stressxx_ana(Ly, F, etas, etap, lambdap, y):
    return 2*lambdap*etap*(vgrad(Ly, F, etas, etap, y))**2

def stressxy_ana(Ly, F, etas, etap, lambdap, y):
    return etap*vgrad(Ly, F, etas, etap, y)


nruns = 4
L = np.full(nruns,30,dtype=int)
F = np.full(nruns,1.48e-5)
lambdap = np.full(nruns,1801.8)
umax = np.empty(nruns)
etas=0.5
etap=0.5

stress = []
stressa = []
vel = []
vela = []
y = []
relerr = np.zeros((nruns,3))

for i in range(nruns):
    L[i] = L[i]*(i+1)
    F[i] = F[i]/(i+1)**3
    lambdap[i] = lambdap[i]*(i+1)**2
    umax[i] = F[i]*L[i]**2/(8.*(etas+etap))
    stress.append( np.zeros((L[i],2)) )
    stressa.append( np.zeros((L[i],2)) )
    vel.append( np.zeros(L[i]) )
    vela.append( np.zeros(L[i]) )
    y.append( np.zeros(L[i]) )

diff = np.zeros((nruns,L[-1]))

for j in range(nruns):
    with open('stress' + str(L[j]) + '.csv','rb') as f:
        reader = csv.reader(f)
        i = 0
        for row in reader:
            stress[j][i,0] = float(np.array(row)[0])
            stress[j][i,1] = float(np.array(row)[1])
            i += 1

    with open('vel' + str(L[j]) + '.csv','rb') as f:
        reader = csv.reader(f)
        i = 0
        for row in reader:
            vel[j][i] = float(np.array(row)[0])
            y[j][i] = float(np.array(row)[6]) - 0.5
            i += 1
            
    vela[j] = parabola(L[j],F[j],etas,etap,y[j])
    stressa[j][:,0] = stressxx_ana(L[j],F[j],etas,etap,lambdap[j],y[j])
    stressa[j][:,1] = stressxy_ana(L[j],F[j],etas,etap,lambdap[j],y[j])

    relerr[j,0] = np.sqrt(np.sum(np.power(vel[j][1:-1] - vela[j][1:-1],2))/np.sum(np.power(vela[j][1:-1],2)))
    relerr[j,1] = np.sqrt(np.sum(np.power(stress[j][1:-1,0] - stressa[j][1:-1,0],2))/np.sum(np.power(stressa[j][1:-1,0],2)))
    relerr[j,2] = np.sqrt(np.sum(np.power(stress[j][1:-1,1] - stressa[j][1:-1,1],2))/np.sum(np.power(stressa[j][1:-1,1],2)))

plotwidth = 5
fig,ax = plt.subplots(figsize = (plotwidth,plotwidth*0.75))

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

ax.plot(L,relerr[:,0],'-o',label=r'$u_x$',color=colors[0],linewidth=1.5)
ax.plot(L,relerr[:,1],'-o',label=r'$\tau_{xx}$',color=colors[1],linewidth=1.5)
ax.plot(L,relerr[:,2],'-o',label=r'$\tau_{xy}$',color=colors[2],linewidth=1.5)
ax.set_xlabel(r'$L_y\ (\delta x)$')
ax.set_ylabel('Relative Error')
ax.legend(loc=0,numpoints=1,frameon=False)
ax.set_xlim(20,140)
ax.set_ylim(0,0.014)
plt.savefig('/home/cstewart/thesis/plots/pois_l2_error_L.pgf', bbox_inches='tight')
