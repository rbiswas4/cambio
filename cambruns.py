#!/usr/bin/env python

import subprocess
import sys
import time

def runcamb(execname , outfile, paramfile,dowait=False):

	starttime = time.time()
	success = False
	of = open(outfile,'a+')
	processstring = execname + ' ' + paramfile
	of.write (processstring + "\n")
	of.close()
	of = open(outfile,'a+')

	child = subprocess.Popen(processstring, 
		shell = True, 
		stdout = subprocess.PIPE, 
		stderr = subprocess.PIPE)
	if (dowait):
		child.wait()
	while True:
		out = child.stdout.readline()	
		of.write(out)
		if out.find("DONE")!=-1:
			success = True
		if child.poll()!=None:
			of.write ("\n")
			break

	of.close()
	endtime = time.time()
	
	elapsed = endtime - starttime
	if success:
		return (0,  elapsed) 
	else:
		return (1 , elapsed) 

#print runcamb('camb', 'out','params.ini')
