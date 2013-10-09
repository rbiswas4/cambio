#!/usr/bin/env python

"""
Set of routines to plot quantities from CAMB output files, as opposed to 
using pycamb to run CAMB 

"""
import matplotlib.pyplot as plt
import numpy as np 

def __loadtransfers(rootname = None, 
	filename = None,
	dirname = None):

	""" 
	returns transfer functions from the CAMB output 

		col	: optional, defaults to 1 (CDM)
			int, col corresponding to quantity of interest
			1 = CDM
			2 = Baryon
			3 = radiation
			4 = neutrino 
			5 = massive neutrino
			6 = Total 
	"""

	if dirname != None:
		rootname = dirname + "/" + rootname
	if filename == None :
		filename  = rootname + "_transfer_out.dat" 

	transfers = np.loadtxt(filename) 


	return transfers

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
		if As == None or ns == None :
			if transfers:
				raise ValueError()
		else:
			PPS =  __PrimordialPS (koverh, ns , As , h ) 	
	else :
		PPS = np.ones(len(koverh))

	
	
	transsq = matterpower / (2.0*np.pi*np.pi* h *h * h  * k * PPS ) 
	trans = np.sqrt(transsq)

	return trans 

def __preparetransfersforcombination ( transfers, 
	collist ):

	"""
	Returns tuples of transfers used in __combinetransfers 
	given a numpy array transfers (as loaded from CAMB for example)  
	and the list of columns collist for transfer functions. 


	args:

	returns:
	"""

	
	koverh  = transfers[:,0] 


	slist = [] 
	for i in range(len(col)) :
		s = np.zeros(shape = (len(koverh) ,2 ))
		s[:,0] = koverh
		s[:,1] = transfers[:,col[1]] 


		slist.append(s)

	return tuple(slist) 


def __combinetransfers ( koverh , transfertuples , f ) :
	""" 
	returns the combined tranfer function of the tuple of transfer
	outputs and their koverh values at the points in the 
	numpy array  koverh. koverh must not extend beyong the common 
	koverh of the transfer outputs. 

	args:
		koverh: arraylike 
			numpy array of k/h values at which the transfer function
			will be returned. The range of this koverh should be 
			common to all the koverh ranges of the transfers 
		transfertuples: tuple of transfer functions, mandatory
			Each element of transfertuples 
			should have a col for koverh, and the transfer output 
			from CAMB (obviously of the same length). Different 
			transfers can have different lengths. 
		f: array like
			array of fractions of each of component
		
	returns:

	"""

		# First find the intersection range of k/h 
		# and store them in minkh and maxkh. 
		# Clip koverh to the range of [minkh, maxkh] 

	minkh = 0. 
	maxkh = 1000.0

	for transfers in transfertuples:

		t = transfers[:,0]
		mint = min(t) 
		maxt = max(t)

		if mint > minkh :
			minkh  = mint 
		if max < maxkh :
			maxkh = maxtra

	kbools = (koverh > minkh) and (koverh < minkh ) 
	koverh = np.array(koverh[kbools])
		
		# Use numpy to interpolate the results
	transferlist = []
	for transfers in transfertuples:
		koverh_native = transfers[:,0]
		transfervals  = transfers[:,1] 
		interpolatedtransfers = np.interp (koverh, koverh_native , 
			transvervals )
		transferlist.append(interpolatedtransfers)

	reqdtransfers = np.array(transferlist) 

	v = reqdtransfers * f 

	return  v 
	
	
def __matterpowerfromtransfersforsinglespecies(
	koverh ,
	transfer , h ,
	As ,
	ns ):

	"""
	args:
		koverh :
		transfer :
		koverh :
		h :
		As :
		ns :
	usage:

	"""
	PPS =  __PrimordialPS (koverh, ns , As , h ) 	

	k = koverh*h 

	matterpower = 2.0*np.pi*np.pi* h *h * h * transfer * transfer * k * PPS  

	res = np.zeros(shape= (len(transfer),2))

	res[:,0] = koverh
	res[:,1] = matterpower

	return res 
	
def __PrimordialPS( koverh , ns , As , h , k0 = 0.05 ):
	"""Returns the primordial power spectrum 

	"""

	k = koverh * h 
	return As* (k/k0 )**(ns -1.)  

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

	if filename == None :
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
	plotpk("test",showtitle=True)	
	plotpk("CM_0.04_0.1175_70.0_-0.725_1.0")
	plotpkresids("CM_0.04_0.1175_70.0_-0.725_1.0","test")

