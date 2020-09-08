import numpy as np
import matplotlib.pyplot as plt
import free_energy3 as FE
from scipy.optimize import fsolve
import time
import ternary
from concurrent import futures

def do_gam(rat):
	return list(FE.tangent_points(rat,xsteps,RZ2,RZ3,0))+list(FE.tangent_points(rat,xsteps,RZ2,RZ3,1))

t0 = time.time()

# We need to specify the charges of the two species (note: the charge ratio is what matters)
# and the range of gamma to search over:
#Z1,Z2,Z3 = 8,26,34
#elem = ('O','Fe','Se')
Z1,Z2,Z3 = 6,8,10
elem = ('C','O','Ne')

Z1,Z2,Z3 = 6,8,26
elem = ('C','O','Fe')

Z1,Z2,Z3 = 6,8,12
elem = ('C','O','Mg')

Z1,Z2,Z3 = 6,8,11
elem = ('C','O','Na')


print(Z1,Z2,Z3,elem)

# charge ratios
RZ2 = 1.0*Z2/Z1
RZ3 = 1.0*Z3/Z1
G1, G2 = 0.8, 1.1*RZ3**(5.0/3.0)

# number of steps in x and gamma
xsteps = 200
gsteps = 100

# scan through gamma and find the tangent points
rat_vec = np.arange(gsteps+1)*(G2-G1)/gsteps + G1

# or specify the gamma ratios you want by hand, e.g.:
#rat_vec = [0.94, 1.05, 1.4, 3.0]

with futures.ProcessPoolExecutor() as executor:
	all_results = executor.map(do_gam, rat_vec)
# single thread:
#all_results = [do_gam(rat) for rat in rat_vec]

print("Time taken = ",time.time()-t0, 'seconds')

# output files (two options here depending on whether Z is integer or not)
fp = open('dat/liquidus_%d_%d_%d.dat' % (Z1,Z2,Z3), 'w')
fp2 = open('dat/solidus_%d_%d_%d.dat' % (Z1,Z2,Z3), 'w')
#fp = open('dat/liquidus_%.3g_%.3g_%.3g.dat' % (Z1,Z2,Z3), 'w')
#fp2 = open('dat/solidus_%.3g_%.3g_%.3g.dat' % (Z1,Z2,Z3), 'w')

for results in all_results:
	results = np.array(results)
	if len(results)>0:
		x1_points = results[:,0]
		x2_points = results[:,1]
		gamma_points = results[:,2]
		point_type = results[:,3]
		x1_dest = results[:,4]
		x2_dest = results[:,5]

		for x1,x2,d1,d2,pt in zip(x1_points,x2_points,x1_dest,x2_dest,point_type):
			if pt < 0:  # liquid point
				print("%lg %lg %lg %lg %lg" % (gamma_points[0], x1, x2, d1, d2), file=fp)
			else:  # solid point
				print("%lg %lg %lg %lg %lg" % (gamma_points[0], x1, x2, d1, d2), file=fp2)

fp.close()
fp2.close()
print("Time taken for output = ",time.time()-t0, 'seconds')
