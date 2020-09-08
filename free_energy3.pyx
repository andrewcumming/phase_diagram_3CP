#cython: boundscheck=True, wraparound=False, cdivision=True
#
# Computes the free energy of a 3 component plasma
#
# python setup3.py build_ext --inplace

cdef extern from "math.h":
	double exp(double m)
	double log(double m)
	double sqrt(double m)

import numpy as np
import matplotlib.pyplot as plt

cdef double gamma_crit = 178.6

cdef double df_OCP(double gamma):
	# the difference between liquid and solid free energies for a 
	# one-component plasma. Implements equation (9) of Medin & Cumming (2010)
	# Note that some significant digits were dropped in the paper, and have
	# been added back here to ensure the free energy is continuous at Gamma=100 and 200
	if gamma<100.0:
		return -0.3683 + 0.00458*(gamma-100.0)
	if gamma>200.0:
		return 0.09384 + 0.00428*(gamma-200.0)
	return -0.003243*gamma + 1.8645*gamma**0.32301 - 1.7748*log(gamma) - 0.2316 + 10.84/gamma

cdef int is_liquid(double gamma1, double x1, double x2, double RZ2, double RZ3):
	if f_liquid(gamma1, x1, x2, RZ2, RZ3) < f_solid(gamma1,x1, x2, RZ2, RZ3):
		return -1
	else:
		return 1
		
cdef double Smix(double x1, double x2, double RZ2, double RZ3):
	cdef double x3, xbar
	# Entropy of mixing
	x3 = 1.0-x1-x2
	if x3<0.0:
		x3 = 0.0
	xbar = x1 + x2*RZ2 + x3*RZ3
	Smix_val = 0.0
	if x1>0.0:
		Smix_val += x1*log(x1/xbar)
	if x2>0.0:
		Smix_val += x2*log(x2*RZ2/xbar)
	if x3>0.0:
		Smix_val += x3*log(x3*RZ3/xbar)
	return Smix_val

cdef double f_min(double gamma1, double x1, double x2, double RZ2, double RZ3):
	cdef double fl, fs, Smix_term
	fl = f_liquid(gamma1, x1, x2, RZ2, RZ3)
	fs = f_solid(gamma1, x1, x2, RZ2, RZ3)
	Smix_term = Smix(x1,x2,RZ2,RZ3)
	if fs<fl:
		return fs + Smix_term
	else:
		return fl + Smix_term
		
def dfdx(double gamma1, double x1, double x2, double RZ2, double RZ3):
	cdef double f1,f2,f3,eps,d1,d2
	eps = 1e-8
	if x1==0.0:
		f1 = f_min(gamma1,x1,x2,RZ2,RZ3)
		f2 = f_min(gamma1,x1+eps,x2,RZ2,RZ3)
		d1 = (f2-f1)/eps
	elif x1==1.0:
		f1 = f_min(gamma1,x1-eps,x2,RZ2,RZ3)
		f2 = f_min(gamma1,x1,x2,RZ2,RZ3)
		d1 = (f2-f1)/eps
	else:
		f1 = f_min(gamma1,x1-eps,x2,RZ2,RZ3)
		f2 = f_min(gamma1,x1+eps,x2,RZ2,RZ3)
		d1 = (f2-f1)/(2.0*eps)
	if x2==0.0:
		f1 = f_min(gamma1,x1,x2,RZ2,RZ3)
		f2 = f_min(gamma1,x1,x2+eps,RZ2,RZ3)
		d2 = (f2-f1)/eps
	elif x2==1.0:
		x1 = x1-eps
		f1 = f_min(gamma1,x1,x2-eps,RZ2,RZ3)
		f2 = f_min(gamma1,x1,x2,RZ2,RZ3)
		d2 = (f2-f1)/eps
	else:
		f1 = f_min(gamma1,x1,x2-eps,RZ2,RZ3)
		f2 = f_min(gamma1,x1,x2+eps,RZ2,RZ3)
		d2 = (f2-f1)/(2.0*eps)
	return (d1,d2)

cdef double f_liquid(double gamma1, double x1, double x2, double RZ2, double RZ3):
	cdef double gamma2, gamma3, x3
	# Returns the liquid free energy (minus the solid linear mixing terms)
	x3 = 1.0-x1-x2
	if x3<0.0:
		x3 = 0.0
	gamma2 = gamma1 * RZ2**(5.0/3.0)
	gamma3 = gamma1 * RZ3**(5.0/3.0)
	# the liquid free energy is linear mixing plus Smix
	return x1*df_OCP(gamma1) + x2*df_OCP(gamma2) + x3*df_OCP(gamma3) #+ Smix(x1,x2,RZ2,RZ3)

cdef double delta_g(double x, double RZ):
	cdef double xx, RZm1, C, denom
	xx = sqrt(x)
	RZm1 = RZ-1.0
	C = 0.05*RZm1**2 / ((1 + 0.64*RZm1)*(1 + 0.5*RZm1**2))	
	denom = 1 + 27.0*RZm1*xx*(xx-0.3)*(xx-0.7)*(xx-1.0)/(1.0+0.1*RZm1)
	return C/denom	

cdef double f_solid(double gamma1, double x1, double x2, double RZ2, double RZ3):
	# Returns the solid free energy (minus the linear mixing terms)
	# First calculate the deviation from linear mixing for the solid
	# following Ogata et al. 1993
	cdef double RZ2m1, RZ3m1, xx, C, denom, x3, gamma2
	x3 = 1.0-x1-x2
	if x3<0.0:
		x3 = 0.0
	gamma2 = gamma1 * RZ2**(5.0/3.0)
	# add the deviation from linear mixing and Smix to get the free energy
	return gamma1*x1*x2*delta_g(x2/(x1+x2),RZ2) + gamma1*x1*x3*delta_g(x3/(x1+x3),RZ3) + gamma2*x2*x3*delta_g(x3/(x2+x3),RZ3/RZ2) #+ Smix(x1,x2,RZ2,RZ3)

def tangent_points(double rat, int xsteps, double RZ2, double RZ3, int direction):
	cdef double gamma1, d2, d1
	cdef int count, lastcount, i, j, n, imax
	# set gamma
	gamma1 = gamma_crit/rat

	# calculate the minimum free energy and derivative as a function of x
	n = xsteps+1
	x1 = (0.0+np.arange(n))/xsteps
	x1vec = np.array([])
	x2vec = np.array([])
	fmin = np.array([])
	df1 = np.array([])
	df2 = np.array([])
	for i in range(n):
		for j in range(n-i):
			if direction == 0:
				i1 = j
				i2 = i
			else:
				i1 = i
				i2 = j
			x1vec = np.append(x1vec,x1[i1])
			x2vec = np.append(x2vec,x1[i2])		
			fmin = np.append(fmin, f_min(gamma1,x1[i1],x1[i2],RZ2,RZ3))
			deriv = dfdx(gamma1,x1[i1],x1[i2],RZ2,RZ3)
			df1 = np.append(df1, deriv[0])
			df2 = np.append(df2, deriv[1])
			#print(i,j,x1[i],x1[j])

	# we loop through different composition values, find the tangent line at each
	# point, and then test how many times the tangent line intersects the Fmin curve
	# When that number drops to zero or increases suddenly from zero, we have 
	# a tangent point
	d1_last = 0.0
	d2_last = 0.0
	for i in range(len(x1vec)):
 		# tangent line  fdiff = flin - fmin
		fdiff = fmin[i]-fmin + df1[i]*(x1vec-x1vec[i]) + df2[i]*(x2vec-x2vec[i])
		# calculate number of intersections
		count = (fdiff>0.0).sum()
		imax = np.argmax(fdiff)
		d1 = x1vec[imax]
		d2 = x2vec[imax]
		#imax = np.argmax(np.delete(fdiff,i))
		#d1 = np.delete(x1vec,i)[imax]
		#d2 = np.delete(x2vec,i)[imax]
		# did this change to or from zero? 
		if (lastcount == 0 and count > 0) or (lastcount > 0 and count == 0):
			# we've found a tangent point
			# try to find the other end of the tie line
			if count == 0:
				d1 = d1_last
				d2 = d2_last
			if (direction == 0 and x1vec[i]>0) or (direction==1 and x2vec[i]>0):
				yield (x1vec[i], x2vec[i], rat, is_liquid(gamma1,x1vec[i],x2vec[i],RZ2,RZ3), d1, d2)			
		lastcount = count
		d1_last = d1
		d2_last = d2

def test():
	# Reproduce the free energy curves shown in Medin & Cumming (2010) Fig.1
	RZ = 34.0/8.0
	gamma1 = gamma_crit/6.0

	x2 = np.arange(99)*0.01 + 0.01
	fL = np.array([f_liquid(gamma1,1.0-x,1e-10,1.0,RZ)+Smix(1.0-x,1e-10,1.0,RZ) for x in x2])
	fS = np.array([f_solid(gamma1,1.0-x,1e-10,1.0,RZ)+Smix(1.0-x,1e-10,1.0,RZ) for x in x2])

	plt.plot(x2,fL)
	plt.plot(x2,fS)
	plt.plot(x2,np.minimum(fS,fL))
	plt.xlabel(r'$x_2$')
	plt.ylabel(r'$\mathrm{Free\ energy}$')
	plt.savefig('free_energy3_test.pdf')
	