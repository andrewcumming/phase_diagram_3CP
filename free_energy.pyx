#cython: boundscheck=True, wraparound=False, cdivision=True

cdef extern from "math.h":
	double exp(double m)
	double log(double m)
	double sqrt(double m)

import numpy as np
import matplotlib.pyplot as plt

_gamma_crit = 178.6
cdef double gamma_crit = 178.6

def _df_OCP(gamma):
	return df_OCP(gamma)

cdef double df_OCP(double gamma):
	# the difference between liquid and solid free energies for a 
	# one-component plasma. Implements equation (9) of Medin & Cumming (2010)
	# Note that some significant digits were dropped in the paper, and have
	# been added back here to ensure the free energy is continuous at Gamma=100 and 200
	#return -0.3683 + 0.0046*(gamma-100.0)
	if gamma<100.0:
		return -0.3683 + 0.00458*(gamma-100.0)
	if gamma>200.0:
		return 0.09384 + 0.00428*(gamma-200.0)
	return -0.003243*gamma + 1.8645*gamma**0.32301 - 1.7748*log(gamma) - 0.2316 + 10.84/gamma

cdef int is_liquid(double gamma1, double x1, double RZ):
	if f_liquid(gamma1,x1,RZ) < f_solid(gamma1,x1,RZ):
		return -1
	else:
		return 1
		
cdef double Smix(double x1,double RZ):
	cdef double x2
	# Entropy of mixing
	x2 = 1.0-x1
	return x1*log(x1/(x1+x2*RZ)) + x2*log(x2*RZ/(x1+x2*RZ))

cdef double f_min(double gamma1, double x1, double RZ):
	cdef double fl, fs
	fl = f_liquid(gamma1, x1, RZ)
	fs = f_solid(gamma1, x1, RZ)
	if fs<fl:
		return fs
	else:
		return fl

cdef double dfdx(double gamma1, double x1, double RZ):
	cdef double f1,f2,eps
	eps = 0.0001
	f1 = f_min(gamma1,x1,RZ)
	f2 = f_min(gamma1,x1+eps,RZ)
	return (f2-f1)/eps

cdef double f_liquid(double gamma1, double x1, double RZ):
	cdef double gamma2
	# Returns the liquid free energy (minus the solid linear mixing terms)
	gamma2 = gamma1 * RZ**(5.0/3.0)
	# the liquid free energy is linear mixing plus Smix
	return x1*df_OCP(gamma1) + (1.0-x1)*df_OCP(gamma2) + Smix(x1,RZ)

def _f_solid(gamma1, x1, RZ):
	return f_solid(gamma1, x1, RZ)

cdef double f_solid(double gamma1, double x1, double RZ):
	# Returns the solid free energy (minus the linear mixing terms)
	# First calculate the deviation from linear mixing for the solid
	# following Ogata et al. 1993
	cdef double RZ1, xx, C, denom
	RZ1 = RZ-1.0
	xx = sqrt(1.0-x1)
	C = 0.05*RZ1**2 / ((1 + 0.64*RZ1)*(1 + 0.5*RZ1**2))	
	denom = 1 + 27.0*RZ1*xx*(xx-0.3)*(xx-0.7)*(xx-1.0)/(1.0+0.1*RZ1)
	#C=0.0
	# add the deviation from linear mixing and Smix to get the free energy
	return gamma1*x1*(1.0-x1)*C/denom + Smix(x1,RZ)

def tangent_points(double rat, int xsteps, double RZ):
	cdef double gamma1
	cdef int count, lastcount, i
	# set gamma
	gamma1 = gamma_crit/rat

	# calculate the minimum free energy and derivative as a function of x
	x1 = (1.0+np.arange(xsteps-1))/xsteps
	fmin = np.array([f_min(gamma1,x,RZ) for x in x1])
	df = np.array([dfdx(gamma1,x,RZ) for x in x1])

	# we loop through different composition values, find the tangent line at each
	# point, and then test how many times the tangent line intersects the Fmin curve
	# When that number drops to zero or increases suddenly from zero, we have 
	# a tangent point
	lastcount = -1
	for i in range(xsteps-1):
		# tangent line
		flin = fmin[i] + df[i]*(x1-x1[i])
		## plot
		##plt.plot(x1,flin)
		##plt.plot(x1,fmin)
		##plt.show()
		# calculate number of intersections
		count = (flin>fmin).sum()
		##print(count)
		# did this change to or from zero? 
		if (lastcount == 0 and count > 0) or (lastcount > 0 and count == 0):
			# we've found a tangent point
			yield (1.0-x1[i], rat, is_liquid(gamma1,x1[i],RZ))
		lastcount = count

def test(gamma=6.0,RZ=4.25):
	# Reproduce the free energy curves shown in Medin & Cumming (2010) Fig.1
	# (note that the figure caption has a mistake, it says that the free energies 
	# have had the liquid linear mixing terms subtracted, but it's actually the solid
	# linear mixing terms)
	#RZ = 34.0/8.0
	#RZ = 8.0/6.0
	#gamma1 = gamma_crit/6.0    # gamma = 6.0 for the Medin & Cumming figure
	gamma1 = gamma_crit/gamma
	x2 = np.arange(499)*0.002 + 0.002
	fL = np.array([f_liquid(gamma1,1.0-x,RZ) for x in x2])
	fS = np.array([f_solid(gamma1,1.0-x,RZ) for x in x2])

	plt.plot(x2,fL)
	plt.plot(x2,fS)
	plt.plot(x2,np.minimum(fS,fL))
	plt.xlabel(r'$x_2$')
	plt.ylabel(r'$\mathrm{Free\ energy}$')
	plt.title(r'$\Gamma = %g$' % (gamma,))
	plt.savefig('free_energy_test.pdf')
