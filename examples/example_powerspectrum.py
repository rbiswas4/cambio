#!/usr/bin/env python

import os
import matplotlib.pyplot as plt
import camb_utils
from camb_utils.cambio import cbpowerspectrum
example_data =  os.path.join(os.path.split(camb_utils.__file__)[0], 'example_data')
print (os.listdir(example_data))
ps = cbpowerspectrum(transferfile=os.path.join(example_data, 'oneh_transfer_out.dat'),
                     Omegacdm=0.3,
                     Omegab=0.05,
                     h=0.71,
                     Omeganu=0.0,
                     As=None,
                     ns=None,
                     koverh=None)

plt.loglog(ps[:, 0], ps[:, 1])
plt.show()
