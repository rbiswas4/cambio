#!/usr/bin/env python

#Show that cbtransfer works as expected with an example dataset
import sys
import os
from camb_utils import cambio
from camb_utils import example_data
import numpy as np
import matplotlib.pyplot as plt

# example_data =  os.path.join(os.path.split(camb_utils.__file__)[0], 'example_data/')
        
fname=os.path.join(example_data, 'oneh_transfer_out.dat')
t = cambio.loadtransfers(filename = fname )
transfertuples = cambio.__preparetransfersforcombination(t,[1,2]) 
#c = cambio.__combinetransfers(transfertuples, f = [0.24 ,0.04])
c = cambio.cbtransfer(fname , 0.24, 0.04)
plt.plot(t[:,0], (0.24* t[:,1] + 0.04 * t[:,2])/0.28, label ="by hand")
plt.plot(c[0], c[1], label =" combine")
print np.shape(c[0])
print np.shape(c[1])
plt.yscale('log')
plt.xscale('log')
plt.legend(loc= "best")
plt.show()
