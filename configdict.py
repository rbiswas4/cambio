#!/usr/bin/env python
def repdictval(fname , paramdict , oname , ignorestrings=['#'], dictdelim = '='):
	"""repdictval (fname , paramdict, oname, ignorestrings=["#], 
	dictdelim ='=') rewrites the file fname to oname but replacing
	all keys in fname according to the dictionary paramdict
	args:
		fname:		mandatory
			filename which is being read in
		oname:		mandatory
			filename to which the output will be written
		paramdict:	mandatory
			dictionary which has the values corresponding 
				to keys in fname
		ignorestrings:	optional, defaults to ['#']
			list of strings, after which the remaining part of
			any line in fname will be ignored in making keys
		dictdelim:	optional, defualts to '-'	
			delimiter used to separate keys, values
			in fname
			
	 
	"""
	f = open(fname, "r")
	line = f.readline()
	i = 0
	w = open(oname, "w")
 

        while line != '': 
                tmp = line.strip()
                if tmp :
                        for st in ignorestrings:

                                tokens = tmp.split(st)
				relevant = tokens[0]
				length =len(tokens)

                        	if len(relevant) >1: 
						#Not a comment line
                                	tp = relevant.split(dictdelim)
                             		key = tp[0].strip()
					val = tp[1].strip()
						#replace val
                                	myval = paramdict[str(key)]
					myline = key + ' ' + dictdelim + ' ' 
					comm =''
					if val != myval:
					 	comm="\n" + ignorestrings[0]+"replaced"
					myline += str(myval) +"\n"+ str(comm) +"\n"  
					w.writelines(myline) 
				else:
						#comment line, so just write 
					#print line
					w.writelines(line)
				
                line=f.readline()
    
        f.close()
	w.close()


	

def builddict(fname,ignorestrings=['#'],dictdelim='='):
	"""builddict (fname) reads in the file with filename
	fname, and builds a dictionary of keys vs values from
	it
	args: 
		fname:		mandatory
			filename from which the dictionary is to be 
				built
		ignorestring: 	optional, defaults to ["#"]
			list of strings, after which the remaining
			part of the line should be ignored. 
		dictdelim:	optional, defaults to '='
			delimiter used to separate keys, values
			in building the dictionary
	returns:
		dictionary of keys and values (in strings)
	usage :
		builddict ( fname)  
	status: 
		Seems to work correctly, tested on CAMB params.ini,
		R. Biswas, July 08, 2012 
	"""
	f = open(fname, "r")
	line = f.readline()
	i = 0
	
	paramdict={}
	while line != '':
		tmp = line.strip()
		if tmp :
			for st in ignorestrings:
				tmp = tmp.split(st)[0]
				if len(tmp) >1:
					tp = tmp.split(dictdelim)
					key = tp[0].strip()
					val = tp[1].strip()
					paramdict[str(key)] = str(val)  
		line=f.readline()
	
	f.close()
	return paramdict

def method4(fname):
	"""method4 (fname) eneric method to tokenize a file and print out 
	a column (of any type) assuming that the file is a 
	consistent table.
	
	args:
	returns:
	example usage:
	status:
	
	"""
	#jfrom cStringIO import StringIO
	#from tokenize import generate_tokens
	import re
	print "Method 4: read in files by line"
	print "and rather than printing out all of it, only print out specific cols "
	f = open(fname,"r")
	line = f.readline()
	i = 0 
	
	while line != '':
		tmp= line.strip()
		if tmp :
			#print tmp
			#tmp = line.strip()
			tmpp = tmp.split()
			#i +=1
			#print len(tmpp)
			if len(tmpp) >1:
				print tmpp[1]
		#tmp = line.split(' ')
		#i += 1
		#tmp = 'sdklsd sdjlks '
		#print len(tmp)
		#if len(tmp) > 1: 
			#print tmp[1]
		line=f.readline()
	
	f.close()
	print "Method 4 done"


