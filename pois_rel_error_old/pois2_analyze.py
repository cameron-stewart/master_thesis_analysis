import csv 
import numpy as np
import matplotlib.pyplot as plt

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
    relerr[j,0] = np.sqrt(np.sum(np.power(vel[j,:] - vela[j,:],2))/np.sum(np.power(vela[j,:],2)))
    relerr[j,1] = np.sqrt(np.sum(np.power(stress[j,:,0] - stressa[j,:,0],2))/np.sum(np.power(stressa[j,:,0],2)))
    relerr[j,2] = np.sqrt(np.sum(np.power(stress[j,:,1] - stressa[j,:,1],2))/np.sum(np.power(stressa[j,:,1],2)))

fig1 = plt.figure(figsize = (8/1.5,6/1.5))
plt.semilogx(wi,relerr[:,0],'o')
plt.xlabel('$\mathrm{Wi}$', size=20)
plt.ylabel('Relative error of $v_x$', size=20)
plt.savefig('/home/cstewart/thesis/figures/pois_vel_l2.eps', format='eps', dpi=1200, bbox_inches='tight')

fig2 = plt.figure(figsize = (8/1.5,6/1.5))
plt.semilogx(wi,relerr[:,1],'o')
plt.xlabel('$\mathrm{Wi}$', size=20)
plt.ylabel(r'Relative error of $\tau_{xx}$', size=20)
plt.axis([1e-3,1e3,3.90e-2,4e-2])
plt.savefig('/home/cstewart/thesis/figures/pois_xx_l2.eps', format='eps', dpi=1200, bbox_inches='tight')

fig3 = plt.figure(figsize = (8/1.5,6/1.5))
plt.semilogx(wi,relerr[:,2],'o')
plt.xlabel('$\mathrm{Wi}$', size=20)
plt.ylabel(r'Relative error of $\tau_{xy}$', size=20)
plt.axis([1e-3,1e3,0.0160,0.0165])
plt.savefig('/home/cstewart/thesis/figures/pois_xy_l2.eps', format='eps', dpi=1200, bbox_inches='tight')
