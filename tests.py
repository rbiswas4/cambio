#!/usr/bin/env python

import camb as cp
import numpy as np
import matplotlib.pyplot as plt

def testlinearityofsigma8withAs():

	"""Creates Plot of sigma8 varying As for z =0, keeping 
	other parameters fixed. Should be linear. 


	"""

		#Setup an array of scalar amplitude values

	As  = np.arange(2.0, 3.0 , 0.1)
	As *= 1.0e-9


		#Initialize sigma8 array 
	Sigma8s = np.zeros(len(As))

	for i in range(len(As)):
		params = {"scalar_amp": As[i]}
		Sigma8s [i]  = cp.sigma8(**params)

	
	plt.plot (As, Sigma8s )
	plt.xlabel(r'$A_s$')
	plt.ylabel(r'$\sigma_8$')
	plt.title('Should be linear')
	plt.tight_layout()
	plt.savefig("sigma8vsAs.pdf")

def testAsforsigma8():

	""" Check that the value of As obtained for a particular 
	sigma8 value indeed returns that value.

	"""

	sig8val = 0.79 
	As =  cp.Asforsigma8(sig8val)
	params = {"scalar_amp": As}
	print As
	print "The next line is the difference between the sigma8 value input \n and the sigma8 value found when the obtained value of As is used. This should \n be small"
	print cp.sigma8(**params) - sig8val
	return 0
	
if __name__=="__main__":

	import camb as cp
	#Should run the next test , commented out for speed.
	#testlinearityofsigma8withAs()
	
