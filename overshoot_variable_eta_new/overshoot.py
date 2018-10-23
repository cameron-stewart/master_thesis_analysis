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
betas = np.array([1,3,5,7,9])
nruns = betas.size

tmax = np.zeros(nruns)
vmax = np.zeros(nruns)
vel = np.zeros((leng,nruns))
x = np.arange(leng)

j=0
for beta in betas:
    with open('vel' + str(beta) + '.dat', 'rb') as f:
        reader = csv.reader(f)
        i=0
        for row in reader:
            for col in row:
                if col != row[-1]:
                    if i < leng:
                        vel[i,j] = float(np.array(row)[i])
                    i+=1
    j +=1
    f.close()

vel = vel/vel[-1,-1]

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

plt.legend(frameon=False)
plt.xlabel("Time Steps" + r'$\ (\delta t)$')
plt.ylabel('Relative Velocity' + r'\ $u_x/u_x^\mathrm{ss}$')
plt.savefig("/home/cstewart/thesis/plots/overshoot_beta.pgf", bbox_inches="tight")

for i in range(nruns):
    vmax[i] = np.amax(vel[:,i])
    tmax[i] = np.argmax(vel[:,i])

plt.figure(figsize = (plotwidth,plotwidth*0.75))
#plt.plot(betas/10.,vmax, 'o')
#plt.axis([0,1,1,1.8])
#plt.plot(betas/10.,tmax, 'o')
#plt.axis([0,1,300,500])
