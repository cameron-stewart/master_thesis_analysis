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


def N(n):
    return (2*n-1)*np.pi

def Alpha_n(N):
    return 1+(1-beta)*El*N*N/4

def Beta_n(alpha_n,N,flag):
    if(flag):
        return np.sqrt(alpha_n*alpha_n - El*N*N)
    else:
        return np.sqrt(El*N*N - alpha_n*alpha_n)

def Gamma_n(N):
    return 1-(1+beta)*El*N*N/4

def G(alpha_n,beta_n,gamma_n,flag,T):
    if(flag):
        return ((1.0 - gamma_n/beta_n)*np.exp(-(alpha_n+beta_n)*T/2)+(1.0 + gamma_n/beta_n)*np.exp((beta_n-alpha_n)*T/2))
    else:
        return 2*np.exp(-alpha_n*T/2)*(np.cos(beta_n*T/2) + (gamma_n/beta_n)*np.sin(beta_n*T/2))

def W(T):
    W_ = 1.5
    for n in range(1,1000):
        N_ = N(n)
        alpha_n = Alpha_n(N_)
        
        if(alpha_n*alpha_n - El*N_*N_ < 0):
            flag_ = False
        else:
            flag_ = True

        beta_n = Beta_n(alpha_n,N_,flag_)
        gamma_n = Gamma_n(N_)
        G_=G(alpha_n,beta_n,gamma_n,flag_,T)

        W_ -= 24*(np.sin(N_/2)/(N_*N_*N_))*G_

    return W_

plotwidth = 6 
leng = 12000
lambdas = np.array([1000,3000,5000,7000,9000])
nruns = lambdas.size

beta = 0.3
C = 0.0000148
eta_0 = 1.0
h = 33.0
rho = 1.0


W_m = h*h*C/(12*eta_0)

Wi = lambdas*W_m/(h/2)
Re = rho*W_m*(h/2)/eta_0
Els = Wi/Re



t1 = np.zeros(nruns)
vmax = np.zeros(nruns)
t2 = np.zeros(nruns)
v2 = np.zeros(nruns)
vel = np.load('numpy.save.npy')
x = np.arange(leng)

vel = vel/(W_m*1.5)
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
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

for i in range(nruns):
    plt.plot(x[:6000:200],vel[:6000:200,i],".",color=colors[i],label=(r'$\lambda_p =\ $' + str(lambdas[i])))

plt.ylim(0.7,1.4)
plt.xlim(0,6000)
plt.legend(frameon=False)
plt.xlabel("Time Steps" + r'$/\delta t$')
plt.ylabel('Relative Velocity' + r'\ $u_x/u_x^\mathrm{ss}$')

relerr = np.zeros(lambdas.size)

for j in range(lambdas.size):
    T = np.linspace(0,12000/lambdas[j],12000)
    El = Els[j]
    W_ = W(T)
    w_ = W_*W_m/(W_m*1.5)
    plt.plot(T*lambdas[j],w_, color=colors[j], linewidth=1.0)
    relerr[j] = np.sqrt(np.sum(np.power(w_ - vel[:,j],2))/np.sum(np.power(w_,2)))

plt.savefig("/home/cstewart/thesis/plots/overshoot_compare_lambda.pgf", bbox_inches="tight")
#plt.savefig("overshoot_compare_lambda.pdf", bbox_inches="tight")


plt.figure(figsize = (plotwidth,plotwidth*0.6))

ax2 = plt.subplot(111)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.get_xaxis().tick_bottom()
ax2.get_yaxis().tick_left()
ax2.set_ylim(0.0,0.1)
plt.plot(lambdas,relerr,marker='.',color=colors[1])
plt.xlabel("Polymer Relaxation Time " + r'$\lambda_p/\delta t$')
plt.ylabel("Relative Error")
plt.savefig("/home/cstewart/thesis/plots/overshoot_error_lambda.pgf", bbox_inches="tight")
#plt.savefig("overshoot_error_lambda.pdf",bbox_inches="tight")

