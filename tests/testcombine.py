#!/usr/bin/env python

#Show that cbtransfer works as expected with an example dataset

import numpy as np
from camb_utils import cambio
import matplotlib.pyplot as plt

location = "../example_data"
fname  = location + "/oneh_transfer_out.dat"
t = cambio.__loadtransfers(filename = fname )
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
