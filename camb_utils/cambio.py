#!/usr/bin/env python

"""
Set of routines to 
	(a) load quantities from CAMB output files
	(b) compute very simple implicit functions in creating the CAMB 
		outputs from one another.
	(c) plot basic quantities from CAMB outputs  
Notes:
	(a) Right now, these are alternative to the pycamb routine. Will worry 
		about integration later. 
	(b) No cosmological calculations should be done in this module 
	(somewhat hard to distinguish from b above). For example, we 
	will not normalize to a sigma8 in this routine, but perform that 
	normalization elsewhere. 

R. Biswas, Fri Dec 13 17:01:22 CST 2013
------------------------------------------------------------------------------

Useful Routines: 
Power Spectrum Calculation:
	- Compute power spectrum of a combination of components in a transfer function.
"""
import matplotlib.pyplot as plt
import numpy as np 




def loadpowerspectrum( powerspectrumfile ):
	
	"""Loads the power spectrum from CAMB into an array with columns 
	of k in units of h / Mpc and power spectrum 

	args:
		powerspectrumfile : string, mandatory
			filename of the power spectrum output from CAMB

	"""
	
	import numpy as np

	pk = np.loadtxt(powerspectrumfile) 

	return pk
def loadtransfers(rootname = None, 
	filename = None,
	dirname = None):

	""" 
	returns an array containing the output of transfer function output
	output file from CAMB 

	args:
		rootname: string,optional , defaults to None
			CAMB convention of rootname
		filename: string, optional, defaults to None
			Note: either rootname or filename need to be 
			provided. filename overrides rootname /dirname
		dirname: string , optional, defaults to None

	returns : numpy array of transfer function output from CAMB. The 
		cols of the array are given by 
			0 = koverh in units of h over Mpc
			1 = CDM
			2 = Baryon
			3 = radiation
			4 = neutrino 
			5 = massive neutrino
			6 = Total 

	"""

	if dirname is not None:
		rootname = dirname + "/" + rootname
	if filename is None :
		filename  = rootname + "_transfer_out.dat" 

	transfers = np.loadtxt(filename) 
	return transfers

def PrimordialPS( koverh , ns , As , h , k0 = 0.05 ,nrun = 0.0):
	"""Returns the primordial power spectrum 

	args:
		koverh :
		ns     : float , mandatory 
			scalar spectral index
		As     : float , mandatory 
		h      : float, mandatory
		k0     : float, optional, defaults to 0.05
			pivot, in units of Mpc^{-1}

		nrun   : float, optional, defaults to 0.0
			

	status:
		Tested withoutrunning. Tests with nrun not done
	notes: 
		Consider moving to cosmodefs

	"""
	#print "PRIMORDIAL ", ns , As , k0 
	k = koverh * h 
	running = (k/k0)**(0.5*np.log(k/k0)*nrun ) 
	return As* (k/k0 )**(ns -1.) * running  
	

def cbtransfer ( transferfile, 
	Omegacdm , 
	Omegab , 
	koverh = None ):

	"""
	Returns the baryon- CDM transfer function at requested koverh values
	from the transfer function file and the values of Omegacdm and Omegab
	provided. If no koverh values are requested, then the transfer function 
	values are returned at the koverh values of the  provided transfer 
	function file.  
	args:
		transferfile: string, mandatory
			absolute path to transfer function file produced by CAMB
		Omegacdm    : float , mandatory
			value of omegacdm used to combine transfer functions

		Omegab      : float mandatory
			value of Omegab used to combine transfer functions
		koverh      : array like ,optional defaults to None
			values of koverh at which the values are requested.
			if None, the transfer functions are returned at the 
			values of koverh in the input transfer function file

	returns: 
		tuple of koverh , Tk of CDM and baryon combined (as in CAMB )

	"""
 	transfers = loadtransfers(rootname = None, 
		filename = transferfile) 
	f = [Omegacdm , Omegab ] 
	tcb =__preparetransfersforcombination(transfers, collist=[1,2])
	tcbcomb = __combinetransfers(tcb , f= f, koverh= koverh)

	
	if koverh is None:
		koverh = transfers[:,0]
	#print "******************"
	#print tcbcomb 
	#print transfers[:,-1]

	return koverh, tcbcomb 

def matterpowerfromtransfersforsinglespecies(
	koverh ,
	transfer , 
	h ,
	As ,
	ns ):

	"""
	returns the power spectrum values corresponding to a set of 
		transfer function values at the values of koverh 
		(comoving k/h in units of h/Mpc)
	args:
		koverh : array like, mandatory (but can be None) 
			values of comoving k/h  in units of h/Mpc as in CAMB 
			cols at which power spectrum values are requested
			
		transfer :  a tuple of (k/h , transfer functions) 
		h :
		As :
		ns :
	returns  : array of shape (numkoverh, 2) , with arr[:,0] = koverh , 
		arr[:,1] = powerspectra
	usage:
		>>> transferout = loadtransfers(filename = "transfer.dat")

	status:

	"""
	#print "khajksd ", As
	#print "As in matterpowerfromtransfersforsinglespecies" , As

	if koverh is None:
		koverh = transfer[0]

	PPS =  PrimordialPS (koverh, ns , As , h ) 	

	k = koverh*h 

	#print type(transfer)
	transferinterp = np.interp(koverh, transfer[0],transfer[1],left = np.nan, right = np.nan)
	#print "shapes" , np.shape(k) , np.shape(koverh), np.shape(transfer), np.shape(transferinterp)
	matterpower = 2.0*np.pi*np.pi* h *h * h * transferinterp * transferinterp * k * PPS  

	res = np.zeros(shape= (len(koverh),2))

	res[:,0] = koverh
	res[:,1] = matterpower

	return res 

def cbpowerspectrum( transferfile, 
	Omegacdm , 
	Omegab , 
	h , 
	Omeganu = 0.0, 
	As = None, 
	ns = None ,
	koverh = None ):
	"""
	Returns the baryon- CDM matter power spectrum using the transfer function
	output, usng the cosmological parameters As, ns, h, Omegab, Omegacdm 

	args:
		As : If As is None, a default value of 2.1e-9 is applied
		     If As >1e-5, then As is assumed to be a ratio to be 
			multiplied to Asdefault to get the correct As
		     If As < 1e-5, then As is assumed to be the real As value


	returns:
		array res, where res[:,0] is koverh and res[:,1] is the power 
		spectrum
	"""
		
	#print "AS in cbpowerspectrum " , As
	
	#Asdefault = 2.1e-9 
	if As is None :
		As = 1.0#Asdefault
	#elif As > 1e-5:
	#	As = Asdefault *As 
	#else:
	#	As = As 
 

	if ns is None :
		ns = 0.963

 	transfers = loadtransfers(rootname = None, 
		filename = transferfile) 

	f = [Omegacdm , Omegab ] 
		#Do as in HACC
	#print "Omeganu ", Omeganu 

		###The lines below are wrong
		###for the calculation follow comments with three ###
		### R. Biswas ,  Mon Mar 31 19:42:09 CDT 2014
	# want rhob/(rhob + rhonu + rhoc )*Tb + rhoc /(rhoc + rhonu +rhob)*Tc
	# = (rhob+ rhoc)/(rhob + rhoc + rhonu) *   (rhob/(rhob + rhoc) + rhoc/(rhoc + rhob))
	# = (1 - fnu) * (f_b*Tb + f_c*c )
		###I want (rhob * Tb + rhoc * Tc )/(rhob + rhoc)
	fnu = Omeganu / (Omeganu + Omegacdm + Omegab )
	fcb = 1.0 - fnu 
	#print f

	tcb =__preparetransfersforcombination(transfers, collist=[1,2])
	
	#f_b Tb + fc Tc 
	tcbcomb = __combinetransfers(tcb , f= f, koverh= koverh)

	#print "*****************"
	#print transfers[:,0]
	#print tcbcomb
	#print "*****************"

	#print "tcbcomb ", type(tcbcomb)
	#print np.shape(tcbcomb) , len(transfers[:,0])
	if koverh ==None:
		koverh = transfers[:,0]
		res = matterpowerfromtransfersforsinglespecies(
		koverh  = koverh,
		transfer = (koverh, tcbcomb), 
		h = h,
		As = As,
		ns  = ns)
	#print "fcb*f", fcb ,f

	#res [:,1] = fcb * fcb*res[:,1]
	return res	
###########################################################################
#########################                 ################################# 
######################### Helper Routines #################################
#########################                 ################################# 
########################################################################### 
def __preparetransfersforcombination ( transfers, 
	collist ):

	"""
	Returns tuples of transfers used in __combinetransfers 
	given a numpy array transfers (as loaded from CAMB for example)  
	and the list of columns collist for transfer functions. 


	args:
		transfers:
			array of transfer functions formed by loading the 
			CAMB transfer function into a numpy array

		collist : list of column numbers

	returns:
		tuples of numpy arrays. Each numpy array has the 0th column 
		koverh while column 1 is the transfer function output from
		CAMB.  
	"""

	
	koverh  = transfers[:,0] 


	slist = [] 
	for i in range(len(collist)) :
		s = np.zeros(shape = (len(koverh) ,2 ))
		s[:,0] = koverh
		s[:,1] = transfers[:,collist[i]] 


		slist.append(s)

	return tuple(slist) 



def __combinetransfers ( transfertuples , f , koverh = None) :
	""" 
	returns the combined tranfer function of the tuple of transfer
	outputs and their koverh values at the points in the 
	numpy array  koverh. koverh must not extend beyond the common 
	koverh of the transfer outputs. 

	args:
		transfertuples: tuple of transfer functions, mandatory
			Each element of transfertuples 
			should have a col for koverh, and the transfer output 
			from CAMB (obviously of the same length). Different 
			transfers can have different lengths. 
		f: array like
			array of fractions of each of component
		koverh: arraylike 
			numpy array of k/h values at which the transfer function
			will be returned. The range of this koverh should be 
			common to all the koverh ranges of the transfers 
		
	returns:
		array of combined transfer function of len equal to that of the
		koverh array

	"""


	koverh_native = transfertuples[0][:,0]

		#Where to interpolate
	if koverh ==None:
		koverh = koverh_native 
	else:
			#If koverh supplied, clip to range of koverh_native
			# ie. we will not extrapolate
		mint = min(koverh_native) 
		maxt = max(koverh_native)
		
		kbools = (koverh > mint) & (koverh < maxt ) 
		koverh = np.array(koverh[kbools])

		#Put (interpolated if k values provided) transfer functions 
		#into an array

	transferlist = []
	for transfers in transfertuples:
		transfervals  = transfers[:,1] 
		interpolatedtransfers = np.interp (koverh, koverh_native , 
			transfervals ,left = np.nan, right= np.nan)
		transferlist.append(interpolatedtransfers)

	reqdtransfers = np.array(transferlist).transpose()

		#Normalize
	f  = np.asarray(f) 
	totfrac = f.sum()	
	fracs = f/totfrac

	v = reqdtransfers*fracs
	#print "all three", np.shape(reqdtransfers), np.shape(fracs), np.shape(v)

	ret =   v.sum(axis=1)
	#print "return shape", np.shape(ret)
	#print "jkhhajkdhask  transfertuples"
	#print transfertuples[0][:,1]
	#print transfertuples[1][:,1]
	#print fracs
	
	ret = transfertuples[0][:,1] *fracs[0]+ transfertuples[1][:,1]*fracs[1]
	#print test
	#ret = test
	return ret


def __matterpowerfromtransfers ( transfers , 
	col , 
	h ,
	As ,
	ns ):

	"""
	return the matter power spectrum from the transfer output of CAMB input 
	as an array . As opposed to plotting from the matter power output of 
	CAMB, this allows one to plot the power spectrum of individual 
	components. 

	usage:
		>>> transfers = cambplots.__loadtransfers(rootname = "m000n0", dirname = "../data/CAMB_outputs/")
		>>> pkfromtransfers = cambplots.__matterpowerfromtransfers(transfers, h = 0.71, As = 2.14e-9, ns = 0.963)
	
	status:
		Seems to match the power spectrum computed directly, but have
		not done a ratio plot requiring interpolation.
		
		>>> pk = np.loadtxt ("../data/CAMB_outputs/m000n0_matterpower.dat")
		>>> plt.loglog (pk[:,0], pk[:,1])
		>>> plt.loglog (pkfromtransfers[:,0], pkfromtransfers[:,1])

	"""
	koverh = transfers[:,0]

	TKtot  = transfers[:,col] 

	res  = __matterpowerfromtransfersforsinglespecies ( koverh , 
		TKtot , 
		h , 
		As , 
		ns )
	return res

def __densitycontrastfrommatterpower(
	koverh ,
	matterpower, 
	h  , 
	As = None ,
	ns = None ,
 	transfers = False ) :

	"""
	return the transfer output of a species from the matter power output
	from the same species. If the parameters ns, As for the primordial
	power spectrum are not supplied, this can calculate the transfer 
	output multiplied by the square root of the primordial power spectrum.

	args:

		koverh:
	
		matterpower:
	
		h:
	
		As :
		ns :
	returns:

	usage:
	status:

	"""
	

	if transfers :
		if As is None or ns is None :
			if transfers:
				raise ValueError()
		else:
			PPS =  PrimordialPS (koverh, ns , As , h ) 	
	else :
		PPS = np.ones(len(koverh))

	
	
	transsq = matterpower / (2.0*np.pi*np.pi* h *h * h  * k * PPS ) 
	trans = np.sqrt(transsq)

	return trans 



	
	
def __getdelta ( transfers ,
	z, 
	omegabh2 , 
	omegach2 , 
	omeganuh2 ,
	H0 ):


	"""

	"""
	h = H0/100.
	Omegab  = omegabh2  / h / h
	Omegac  = omegach2  / h / h 
	Omegan  = omeganuh2 / h / h 


	rhob =  Omegab* (1. + z )**3.0
	rhocdm = Omegac* (1. + z)**3.0  
	rhonm = Omegan * (1. + z) **3.0 

	return 0


def plotpk(rootname,
	filename = None, 
	dirname = '',
	color = 'Black', 
	linestyle ='-',
	labels = True, 
	legs = "", 
	title = ""):

	"""Plots matter power spectrum from CAMB output "root_matterpower.dat"
	assuming single redshift output.

	args: 
		rootname: mandatory
			string, output_root of CAMB params.ini file
			eg. "test"
		filename: optional, string 
			overrides the rootname to get the filename of the 
			matter power spectrum rather than use the rootname. 
		dirname : optional, defaults to curren directory
			string, directory where the output of CAMB is stored
			eg. "/home/rbiswas/doc/camboutput/"
		color	: optional, defaults to "Black"
			string, color of plot 
			eg. 'Black'
		linestyle: optional, defaults to '-'
			string, linestyle of plot
			eg. '-'
		labels	: optional, defaults to True
			Bool, Show x, y labels
			eg. False
		legs	: optional, defualts to "", no label for legends
			string, value for label to be used in legend
			eg. "version X"
		title	: optional, defaults to "", no title
			string, value for title of plot
	returns	: 
		0, if successful
		exceptions not defined
	example usage	: 

	status	:seems to work as advertised,
		R.Biswas, July 14, 2012
	"""
	if dirname != '':
		rootname = dirname + "/" + rootname

	if filename is None :
		filename = rootname + "_matterpower.dat"
	data = np.loadtxt(filename)

	if legs !='':
		plt.plot (data[:,0],data[:,1],linestyle,color =color,label =legs)
	else:
		plt.plot (data[:,0],data[:,1],linestyle,color =color)
		
	plt.xscale("log")
	plt.yscale("log")

	if labels==True:
		plt.xlabel("k/h Mpc^{-1}")
		plt.ylabel("P(k)")

	if title!="":
		plt.title(title)
	
	return 0	
	
def plottk(rootname,
	dirname = '',
	col = 1, 
	color = 'Black', 
	linestyle ='-',
	labels = True, 
	legs = "", 
	title = ""):

	"""Plots transfer functions from CAMB output "root_transfer.dat"
	assuming single redshift output.

	args: 
		rootname: mandatory
			string, output_root of CAMB params.ini file
			eg. "test"
		dirname : optional, defaults to curren directory
			string, directory where the output of CAMB is stored
			eg. "/home/rbiswas/doc/camboutput/"
		col	: optional, defaults to 1 (CDM)
			int, col corresponding to quantity of interest
			1 = CDM
			2 = Baryon
			3 = radiation
			4 = neutrino 
			5 = massive neutrino
			6 = Total 
		color	: optional, defaults to "Black"
			string, color of plot 
			eg. 'Black'
		linestyle: optional, defaults to '-'
			string, linestyle of plot
			eg. '-'
		labels	: optional, defaults to True
			Bool, Show x, y labels
			eg. False
		legs	: optional, defualts to "", no label for legends
			string, value for label to be used in legend
			eg. "version X"
		title	: optional, defaults to "", no title
			string, value for title of plot
	returns	: 
		0, if successful
		exceptions not defined
	example usage	: 

	status	:seems to work as advertised,
		R.Biswas, July 14, 2012
	"""
	if dirname != '':
		rootname = dirname + "/" + rootname

	data = np.loadtxt(rootname + "_transfer_out.dat")

	if legs !='':
		plt.plot (data[:,0],data[:,col],linestyle,color =color,label =legs)
	else:
		plt.plot (data[:,0],data[:,col],linestyle,color =color)
		
	plt.xscale("log")
	#plt.yscale("log")

	if labels==True:
		plt.xlabel("k/h Mpc^{-1}")
		plt.ylabel("T(k)")

	if title!="":
		plt.title(title)
	
	return 0	
def plotpkresids(root_fid, 
	root_test, 
	dirname = '' ,
	color = 'Black', 
	linestyle ='-',
	labels=True, 
	epsilon = 0.001, 
	title = '',legends=''):
	"""Plots fractional residuals of matter power spectrum from two 
	CAMB outputs (root_test - root_fid)/root_fid assuming single 
	redshift output.

	args: 
		root_fid: mandatory
			string, output_root of fiducial CAMB params.ini file
			eg. "test1"
		root_test:mandatory
			string, output_root of test CAMB params.ini file
			eg. "test2"
		dirname : optional, defaults to curren directory
			string, directory where the output of CAMB is stored
			eg. "/home/rbiswas/doc/camboutput/"
		color	: optional, defaults to "Black"
			string, color of plot 
			eg. 'Black'
		linestyle: optional, defaults to '-'
			string, linestyle of plot
			eg. '-'
		labels	: optional, defaults to True
			Bool, Show x, y labels
			eg. False
		legends	: optional, defualts to "", no label for legends
			string, value for label to be used in legend
			eg. "version X"
		epsilon	: optional, defaults to 0.001
			checks that the test and fiducial files containing
			the matter power spectrum have the same values of k/h	
			to within epsilon. If this is untrue, the plot routines
			fails
		title	: optional, defaults to "", no title
			string, value for title of plot
	returns	: 
		0, if successful
		1, if the values of k/h in the test and fiducial matter 
			power spectrum files do not match up to within
			epsilon  = 0.001
	example usage	: 

	status	:seems to work as advertised,
		R.Biswas, July 14, 2012
	"""

	if dirname != '':
		root_test = dirname + "/" + root_test
		root_fid = dirname + "/" + root_fid

	datafid = np.loadtxt(root_fid + "_matterpower.dat")
	datatest = np.loadtxt(root_test + "_matterpower.dat")

	if len(datafid[:,0]) !=len(datatest[:,0]):
		return 1
	if any(abs(datatest[:,0] -datafid[:,0])) >epsilon:
		return 1
	if legends !='':
		plt.plot(datafid[:,0],
			(datatest[:,1]-datafid[:,1])/datafid[:,1],
			ls = linestyle,
			color = color,
			label = legends)
	else:
		plt.plot(datafid[:,0],
			(datatest[:,1]-datafid[:,1])/datafid[:,1],
			color = color ,
			ls = linestyle)
	plt.xscale("log")

	plt.axhline( color = 'Black')
	if labels==True:
		plt.xlabel("k/h Mpc^{-1}")
		plt.ylabel("$\Delta P(k)/P(k)$")
	if title != '':
		plt.title(title)	

	return 0
def crossingz( w0 , wa ):
	"""Returns the redshift at which w  crosses -1
	given the w0 and wa for the CPL parametrization	
	w = w0 + (1 - a)wa

	args: 
		w0: mandatory
			float, w0 value
		wa: mandatory
			float, wa value
	returns:
		z:	
			float, redshift at which crossing happens
			Negative z implies no crossing
	status:
	"""	
	a = (w0 + wa +1.)/wa
	z = 1.0/a - 1.0
	return z

if __name__=="__main__":

		
	import matplotlib.pyplot as plt
	ps = cbpowerspectrum(
		transferfile = "example_data/oneh_transfer_out.dat",
		Omegacdm = 0.3, 
		Omegab = 0.05, 
		h = 0.71, 
		Omeganu = 0.0, 
		As = None, 
		ns = None ,
		koverh = None )

	#print type(ps)
	plt.loglog ( ps[:,0], ps[:,1])
	plt.show()
	
	
	#plotpk("test",showtitle=True)	
	#plotpk("CM_0.04_0.1175_70.0_-0.725_1.0")
	#plotpkresids("CM_0.04_0.1175_70.0_-0.725_1.0","test")

