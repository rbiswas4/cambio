#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np 
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
		filename = rootname + "+matterpower.dat")
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
	a = (w0 + wa +1)/wa
	z = 1.0/a - 1.0
	return z

if __name__=="__main__":
	plotpk("test",showtitle=True)	
	plotpk("CM_0.04_0.1175_70.0_-0.725_1.0")
	plotpkresids("CM_0.04_0.1175_70.0_-0.725_1.0","test")

