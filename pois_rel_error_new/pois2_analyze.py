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

def stressxx(Ly, F, etas, etap, lambdap, y):
    return 2*lambdap*etap*(vgrad(Ly, F, etas, etap, y))**2

def stressxy(Ly, F, etas, etap, lambdap, y):
    return etap*vgrad(Ly, F, etas, etap, y)

lambdap = np.array([46,215,1000,4642,21544,100000,464170,2154440,10000000])
L = 30
etap = 0.5
etas=0.5
F = 1e-5
nruns = lambdap.size

stress = np.zeros((nruns,L,3))
stressa = np.zeros((nruns,L,3))
vel = np.zeros((nruns,L))
vela = np.zeros((nruns,L))
wi = np.zeros((nruns))
y = np.zeros((L))
relerr = np.zeros((nruns,4))
diff = np.zeros((nruns,L))

for j in range(nruns):
    with open('stress_' + str(lambdap[j]) + '.csv','rb') as f:
        reader = csv.reader(f)
        i = 0
        for row in reader:
            stress[j,i,0] = float(np.array(row)[0])
            stress[j,i,1] = float(np.array(row)[1])
            stress[j,i,2] = float(np.array(row)[4])
            i += 1

    with open('velocity_' + str(lambdap[j]) + '.csv','rb') as f:
        reader = csv.reader(f)
        i = 0
        for row in reader:
            vel[j,i] = float(np.array(row)[0])
            if(j==0):
                y[i] = float(np.array(row)[6]) - 0.5
            vela[j,i] = parabola(L,F,etas,etap,y[i])
            stressa[j,i,0] = stressxx(L,F,etas,etap,lambdap[j],y[i])
            stressa[j,i,1] = stressxy(L,F,etas,etap,lambdap[j],y[i])
            stressa[j,i,2] = 0.0

            i += 1
            
    wi[j] = lambdap[j]*(parabola(L,F,etas,etap, L/2))/L
    relerr[j,0] = np.sqrt(np.sum(np.power(vel[j,1:-1] - vela[j,1:-1],2))/np.sum(np.power(vela[j,1:-1],2)))
    relerr[j,1] = np.sqrt(np.sum(np.power(stress[j,1:-1,0] - stressa[j,1:-1,0],2))/np.sum(np.power(stressa[j,1:-1,0],2)))
    relerr[j,2] = np.sqrt(np.sum(np.power(stress[j,1:-1,1] - stressa[j,1:-1,1],2))/np.sum(np.power(stressa[j,1:-1,1],2)))

fig1 = plt.figure(figsize = (plotwidth,plotwidth*0.75))

ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

plt.semilogx(wi,relerr[:,0],'o',label=r'$u_x$',color=colors[0])
plt.xlabel('Weissenberg Number' + '$\ \mathrm{Wi}$')
plt.ylabel('Relative Error')
plt.semilogx(wi,relerr[:,1],'o',label=r'$\tau_{xx}$',color=colors[1])
plt.semilogx(wi,relerr[:,2],'o',label=r'$\tau_{xy}$',color=colors[2])
plt.legend(loc=0,numpoints=1)
plt.axis([1e-3,1e3,0.0,0.014])
plt.savefig('/home/cstewart/thesis/plots/pois_l2_error_wi.pgf', bbox_inches='tight')
#plt.savefig('pois_l2_error_wi.pdf', bbox_inches='tight')




