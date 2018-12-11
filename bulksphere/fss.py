import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt;

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

plotwidth = 5

def line(x, m, b):
    return (m*x + b) 

#x = np.array([1./60., 1/66., 1./78., 1./90., 1./108., 1./120., 1./132.])
#y = np.array([0.026074, 0.0267465, 0.027874, 0.0289303, 0.0301774, 0.0308309, 0.0313743])
x = np.array([1./78., 1./90., 1./108., 1./120., 1./132.])
y = np.array([0.027874, 0.0289303, 0.0301774, 0.0308309, 0.0313743])
vpredict = 0.036173
xfit = np.linspace(0,0.02,50)

popt,pcov = opt.curve_fit(line,x,y)

meas_str = "%.3f" % (popt[1]/vpredict)
rel_diff = np.abs((popt[1]-vpredict)/vpredict)
dif_str = "%.1f" % (rel_diff*100.)

yfit  = line(xfit, popt[0], popt[1])
yfit1 = line(xfit, popt[0] + pcov[0,0]**0.5, popt[1] + pcov[1,1]**0.5)
yfit2 = line(xfit, popt[0] - pcov[0,0]**0.5, popt[1] - pcov[1,1]**0.5)

fig, ax = plt.subplots(figsize = (plotwidth,plotwidth*0.75))
ax.plot(x,y/vpredict,'o',label='Measured',color=colors[0])
ax.plot(xfit,yfit/vpredict,color=colors[3], label='Fit')
ax.set_xlabel(r'$1/L/(1/\delta x)$')
ax.set_ylabel('Relative Velocity ' + r'$\ v_t/v_\mathrm{St}\ $')
ax.legend(numpoints=1, frameon=False)
ax.set_xlim(0.0, 0.02)
ax.yaxis.tick_left()
ax.xaxis.tick_bottom()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.text(0.002, 0.028/vpredict,r'$v_t^\mathrm{fss}/v_{St} = $' + meas_str)
plt.text(0.002, 0.0265/vpredict,'Percent Diff.' + r'$\ =\ $' + dif_str +r'$\%$')
plt.text(0.0022, 0.997,r'$v_t^\mathrm{fss}$')
plt.xticks([0.0,0.005,0.01,0.015,0.02])
ax.arrow( 0.002, 1.006, -0.0012, 0., fc="k", ec="k", head_width=0.008, head_length=0.0005, width = 0.0001)
plt.fill_between(xfit, yfit1/vpredict, yfit2/vpredict, facecolor=colors[2], alpha=0.9)
plt.savefig('/home/cstewart/thesis/plots/sphere_bulk_fss.pgf', bbox_inches='tight')
#plt.savefig('sphere_bulk_fss.pdf', bbox_inches='tight')
