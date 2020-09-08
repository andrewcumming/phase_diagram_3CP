import numpy as np
import matplotlib.pyplot as plt
import free_energy as FE
from scipy.optimize import fsolve
import time
from concurrent import futures

def do_gam(rat):
	return list(FE.tangent_points(rat,xsteps,RZ))

def calculate_x_solid():
	# compute the solid composition that freezes out at each liquid composition
	# x_liq, gam_liq, x_sol, gam_sol are already defined
	x_solid = np.array([])
	for i,xl in enumerate(x_liq):
		gam = gam_liq[i]
		xsvals = np.sort(x_sol[gam_sol == gam])
		if len(xsvals)>1:   # if there is more than one option for the solid
							# composition, we need to decide which one to take
			# is gamma increasing or decreasing with x?
			if xl<=x_liq[1]:
				dgam = gam_liq[2]-gam_liq[0]
			else:
				dgam = gam-gam_liq[i-2]
			if dgam>0:  # gamma is increasing with x, take the one to the right
				if xsvals[0]<xl:
					x_to_add = xsvals[1]
				else:
					x_to_add = xsvals[0]
			else:  # gamma is decreasing with x, take the one to the left
				x_to_add = xsvals[0]
		elif len(xsvals)<1:   # if there is no solid composition, use the liquid
			x_to_add = xl
		else:  # if there is only one solid composition, use that one
			x_to_add = xsvals[0]
		x_solid = np.append(x_solid,x_to_add)
	if x_solid[0]>0.0:  # if we get to the left of the diagram, freeze out the light solid
	#	x_solid[x_liq<=0.03]=0.0;
		x_solid[0]=0.0;
	return x_solid


t0 = time.time()

# We need to specify the charges of the two species (note: the charge ratio is what matters)
# and the range of gamma to search over:
#Z1,Z2,G1,G2 = 26,34,0.9,1.7
#Z1,Z2,G1,G2 = 8,34,0.1,12.0
#Z1,Z2,G1,G2 = 2,6,0.1,12.0
#Z1,Z2,G1,G2 = 8,26,0.7,12.0
#Z1,Z2,G1,G2 = 3,4,0.1,1.7
#Z1,Z2,G1,G2 = 2,3,0.8,2.0
#Z1,Z2,G1,G2 =  12,36,0.7,6.4
#Z1,Z2,G1,G2 = 20,40,0.1,3.3
#Z1,Z2,G1,G2 = 3,13,0.2,12.0
#Z1,Z2,G1,G2 = 6,10,0.7,3.0  # CNe
Z1,Z2,G1,G2 = 6,8,0.9,1.7    # CO
#Z1,Z2,G1,G2 = 8,10,0.9,1.5   # ONe
#Z1,Z2,G1,G2 = 3,5,0.7,2.4

# charge ratio
RZ = 1.0*Z2/Z1

# number of steps in x and gamma
xsteps = 500
gsteps = 500

# scan through gamma and find the tangent points
rat_vec = np.arange(gsteps+1)*(G2-G1)/gsteps + G1
#with futures.ProcessPoolExecutor() as executor:
#	results = executor.map(do_gam, rat_vec)
#single process use:
results = [do_gam(rat) for rat in rat_vec]

# extract results
results = np.array([x for res in results for x in res])
x_points = results[:,0]
gamma_points = results[:,1]
point_type = results[:,2]

# sort the points in order of x
if 1:
	ind = np.argsort(x_points)
	x_points = x_points[ind]
	gamma_points = gamma_points[ind]
	point_type = point_type[ind]

# extract the liquidus and solidus curves
x_liq = x_points[point_type<0]
gam_liq = gamma_points[point_type<0]
x_sol = x_points[point_type>0]
gam_sol = gamma_points[point_type>0]

# compute the solid composition that freezes out at each liquid composition
x_solid = calculate_x_solid()

# check the execution time
print("Time taken=",time.time()-t0, 's')

# output to files
np.savetxt('dat/liquidus_%d_%d.dat' % (Z2,Z1), np.c_[x_liq,gam_liq,x_solid])
np.savetxt('dat/solidus_%d_%d.dat' % (Z2,Z1), np.c_[x_sol,gam_sol])

# plot phase diagram
plt.scatter(x_sol,gam_sol,s=4)
plt.scatter(x_liq,gam_liq,s=4)
plt.xlim((0.0,1.0))
plt.xlabel(r'$x_2$')
plt.ylabel(r'$\Gamma_{\rm crit}/\Gamma_1$')

#plt.annotate(r'$\mathrm{O}/\mathrm{Ne}$',size=16,xy=(0.2,0.9*(max(gamma_points)-min(gamma_points))+min(gamma_points)))
plt.annotate(r'$R_Z=%d/%d$' % (Z2,Z1),xy=(0.2,0.9*(max(gamma_points)-min(gamma_points))+min(gamma_points)))

plt.savefig('phase_diagram_%d_%d.pdf' % (Z2,Z1))

