#!/usr/bin/env python

"""
Module of functions to obtain quantities from CAMB easily using the pycamb 
class. 
R. Biswas,  Sun Aug 25 19:27:59 CDT 2013

Requirements: pycamb 

"""

import numpy as np 
import pycamb as pyc
import scipy.interpolate as sip 
from astropy.cosmology import Planck13 as cosmo

def PK ( koverh , redshiftlist = [0] ,
	maxk = 1.0 ,
	logk_spacing =0.02 ,
	**params ) :


	kh, PK , sigma8 = pyc.matter_power (redshifts = redshiftlist, 
		maxk = maxk , 
		logk_spacing = logk_spacing ,  
		get_sigma8 = True, 
		**params )
	return  np.interp( koverh, kh, PK )

def sigma_r (
	R = 8.0 ,
	redshiftlist = [0] , 
	kmin = 1.0e-10 , 
	maxk = 1.0 , 
	logk_spacing = 0.02 ,
	tol = 4.8e-4 ,
	rtol = 4.8e-4, 
	**params ):


	import scipy.integrate as si 


	integrand = lambda k : sigma_r_integrand (k, R, redshiftlist = redshiftlist , maxk = maxk , logk_spacing = logk_spacing , **params) 

	sigmasq = si.romberg ( integrand , a = kmin ,b = maxk , tol=tol, 
		rtol = rtol , show=False, divmax=10, vec_func=False)

	return np.sqrt(sigmasq)



def sigma_r_integrand(k , 
	R  = 8.0 ,
	redshiftlist = [0], 
	maxk = 1.0 ,
	logk_spacing = 0.02 , 
	**params ):


	""" Obtain the sigma(R) values for the redshift list provided 

	"""

	h  = cosmo.H0/100.0  
	kh = k/h 
	x = kh* R

	W = 3*(np.sin(x) - (x)*np.cos(x))/(x)**3

	integrand = k * k * W *W  / 2.0 /np.pi/ np.pi
	integrand *= PK ( kh , redshiftlist = redshiftlist,  
		maxk = maxk , logk_spacing =logk_spacing ,  **params ) 

	return integrand	


def sigma8(redshiftlist = [0] , 
	maxk = 1.0 , 
	logk_spacing = 0.02, 
	**params ):

	"""
	returns the value of sigma8 at the redshift specified 
	for the set of parameters passed in. 

	args:
		redshiftlist: list of floats, optional, defaults to [0]
		maxk : float, optional, defaults to 1.0 
			max value of k to which the matter power spectrum
			is calculated
		logk_spacing: float, optional , defaults to 0.02
			spacing in log(k *Mpc/h) 
		**params: dictionary/keyword argument for CAMB parameters

	returns:
		float, sigma8 for the parameters passed 
	
	example usage:
		>>> import camb as cp
		>>> params = {"scalar_amp":2.3e-9}
		>>> cp.sigma8(**params)
		>>> 0.8339590551185299
		>>> #This value may change if the defaults of pycamb change
		

	status:
		Tested to work correctly for single redshifts. 
		See testlinearityofsigma8withAs() in tests.

		R. Biswas, Mon Aug 26 10:31:00 CDT 2013

	"""
	sigma8  = pyc.matter_power(redshifts = redshiftlist, 
		maxk = maxk , 
		logk_spacing = logk_spacing ,  
		get_sigma8 = True, 
		**params )[-1]

	return float(sigma8)
		

def __getsigma8(As , params):

	params["scalar_amp"] = As
	
	return sigma8( **params)  

def Asforsigma8(sigma8val , Asmin = 1.1e-9, Asmax = 5.1e-9, **params):

	"""Returns the value of As, the amplitude of the primordial scalar
	fluctuations that produces a sigma8 value of sigma8val when all other
	cosmological parameter values are set to those in **params. If **params 
	is None , the parameters are fixed to the default values.

	args:
		sigma8val:
			mandatory, float
			Value of sigma8 for which As is wanted

		Asmin :
			optional, float , defaults to 1.1e-9
			Minimum of range of As searched by bisection 
			method

		Asmax :
			optional, float , defaults to 5.1e-9
			Maximum of range of As searched by bisection 
			method

		params :
			Dictionary of parameter values

	returns:
		float, As value for the set of parameters defined in params,
			+ default parameters set in CAMB setdefaults

	example usage:
		>>> sig8val = 0.79 
		>>> As =  cp.Asforsigma8(sig8val)
		>>> print As 
		>>> 2.0638671875e-09
		>>> # The values may change if the default parameters change


	status:
		Tested, and seems to work giving a difference of order 10^{-6}
		between the input sigma8 and the sigma8 for the As obtained. 
		See testAsforsigma8() in tests

		R. Biswas,  Mon Aug 26 10:28:50 CDT 2013

	"""
	import scipy.optimize as opt


	# First get the sigma8 value for the passed parameters. 
	#params["scalar_amp"] = Asguess

	constr = lambda As: __getsigma8( As , params ) -sigma8val
	sig8  = opt.bisect(constr, Asmin , Asmax )
	return sig8

if __name__ == "__main__":

	import numpy as np
	import matplotlib.pyplot as plt
	k = np.arange(1.0e-4, 1.0, 0.001)
	plt.loglog (k , PK(k))
	
	plt.figure()
	plt.loglog (k ,sigma_r_integrand(k))
	plt.ylim(ymax = 0.5)
	plt.show()
	print sigma_r(R = 8.0*0.7)

	#print sigma_r (R = 8) - sigma8()
