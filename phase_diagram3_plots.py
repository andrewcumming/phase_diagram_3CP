import numpy as np
import matplotlib.pyplot as plt
import free_energy3 as FE
from scipy.optimize import fsolve
import time
import ternary
from concurrent import futures
#import matplotlib
#del matplotlib.font_manager.weight_dict['roman']
#atplotlib.font_manager._rebuild()
#    matplotlib.rcParams['font.serif'] = "Times New Roman"
    # Then, "ALWAYS use sans-serif fonts"
#    matplotlib.rcParams['font.family'] = "serif"

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Times New Roman']})
rc('text', usetex=True)

t0 = time.time()

# We need to specify the charges of the two species (note: the charge ratio is what matters)
# and the range of gamma to search over:
#Z1,Z2,Z3 = 8,26,34
#elem = ('O','Fe','Se')

Z1,Z2,Z3 = 6,8,10
elem = ('C','O','Ne')

#Z1,Z2,Z3 = 6,8,26
#elem = ('C','O','Fe')

#Z1,Z2,Z3 = 6,8,12
#elem = ('C','O','Mg')

#Z1,Z2,Z3 = 6,8,11
#elem = ('C','O','Na')

print(Z1,Z2,Z3,elem)

# charge ratios
RZ2 = 1.0*Z2/Z1
RZ3 = 1.0*Z3/Z1

# read in results
liquidus = np.loadtxt('dat/liquidus_%d_%d_%d.dat' % (Z1,Z2,Z3))
solidus = np.loadtxt('dat/solidus_%d_%d_%d.dat' % (Z1,Z2,Z3))
#liquidus = np.loadtxt('dat/liquidus_%.3g_%.3g_%.3g.dat' % (Z1,Z2,Z3))
#solidus = np.loadtxt('dat/solidus_%.3g_%.3g_%.3g.dat' % (Z1,Z2,Z3))

# extract the gamma values
gam_vec = np.array([])
for gam in liquidus[:,0]:
    if len(gam_vec)==0:
        gam_vec = np.append(gam_vec,gam)
    elif gam != gam_vec[-1]:
        gam_vec = np.append(gam_vec,gam)

count = 0

# make a figure for each gamma
for gam in gam_vec:
    ind = liquidus[:,0] == gam
    points_liq = [(x,y,z) for x,y,z in zip(liquidus[ind,1],liquidus[ind,2],1.0-liquidus[ind,1]-liquidus[ind,2])]
    dest_liq = [(x,y,z) for x,y,z in zip(liquidus[ind,3],liquidus[ind,4],1.0-liquidus[ind,3]-liquidus[ind,4])]
    ind = solidus[:,0] == gam
    points_sol = [(x,y,z) for x,y,z in zip(solidus[ind,1],solidus[ind,2],1.0-solidus[ind,1]-solidus[ind,2])]
    dest_sol = [(x,y,z) for x,y,z in zip(solidus[ind,3],solidus[ind,4],1.0-solidus[ind,3]-solidus[ind,4])]
    
    figure, tax = ternary.figure()
    ax = tax.get_axes()
    ax.axis('off')
    
    tax.boundary(linewidth=0.5)
    tax.gridlines(multiple=0.25, color='black', linewidth=0.5, zorder=0)#, color="blue")

    plt.annotate(r'$Z_1=%d$' % (Z1,), fontsize=14, xy=(-0.03,0.85), fontname='Times New Roman')
    plt.annotate(r'$Z_2=%d$' % (Z2,), fontsize=14, xy=(-0.03,0.78), fontname='Times New Roman')
    plt.annotate(r'$Z_3=%d$' % (Z3,), fontsize=14, xy=(-0.03,0.71), fontname='Times New Roman')
    plt.annotate(r'$\Gamma_1 = {:.1f}$'.format(178.6/gam), fontsize=14, xy=(0.78,0.75),
            fontname='serif')

    if len(points_liq)>0:
        for x,d in zip(points_liq,dest_liq):
            tax.line(x, d, linewidth=0.5, marker='',  linestyle="-", color='C2')    
        tax.scatter(points_liq,s=6,c='C1')#, marker='s', color='blue')#, label="Red Squares")

    if len(points_sol)>0:
        tax.scatter(points_sol,s=6, c='C0')#, marker='s', color='red')#, label="Red Squares")
        #for x,d in zip(points_sol,dest_sol):
        #    tax.line(x, d, linewidth=0.5, marker='', color='red', linestyle="-")

    x1 = np.arange(101)*0.01
    x2 = ( gam - RZ3**(5.0/3.0) + (RZ3**(5.0/3.0)-1.0)*x1 ) / (RZ2**(5.0/3.0)-RZ3**(5.0/3.0))
    ind = (x2>=0.0)
    x1 = x1[ind]
    x2 = x2[ind]
    ind = (x2<=1.0-x1)
    x1 = x1[ind]
    x2 = x2[ind]
    x3 = 1.0-x2-x1
    #if len(x1)>0:
    #    tax.plot(zip(x1,x2,x3), linewidth=1,linestyle=':', c='k')
    
    x2 = np.arange(101)*0.01
    x3 = 0.05 + 0.1*x2
    ind = (x2>=0.0)
    x3 = x3[ind]
    x2 = x2[ind]
    ind = (x2<=1.0-x3)
    x3 = x3[ind]
    x2 = x2[ind]
    x1 = 1.0-x2-x3
    #if len(x1)>0:
    #    tax.plot(zip(x1,x2,x3), linewidth=1,linestyle='--')
    
    #tax.left_axis_label(r'$Z_3/Z_1=$'+elem[2], fontsize=10)
    #tax.right_axis_label(r'$Z_2/Z_1=$'+elem[1], fontsize=10)
    #tax.bottom_axis_label(r'$Z_1$', fontsize=10)
    tax.left_axis_label(r'$x_3$', fontsize=16, offset=0.15, fontname='Times New Roman')
    tax.right_axis_label(r'$x_2$', fontsize=16, offset=0.18, fontname='Times New Roman')
    tax.bottom_axis_label(r'$x_1$', fontsize=16, offset=0.12, fontname='Times New Roman')
    tick_formats = {'b': "%.2f", 'r': "%.2f", 'l': "%.2f"}
    ax = tax.get_axes()
    ax.tick_params(direction='out', pad=-50)

    tax.ticks(axis='lbr', linewidth=0.5, multiple=0.25, tick_formats=tick_formats, fontsize=12, offset=0.03) #, fontsize=6)
    
    tax._redraw_labels()
     
    count +=1
    #plt.savefig('out3/phase_diagram_%04d.png' % (count,), bbox_inches='tight', dpi=300)
    plt.savefig('out3_C_O_Ne/phase_diagram_%04d.pdf' % (count,), bbox_inches='tight')
    
    tax.close()