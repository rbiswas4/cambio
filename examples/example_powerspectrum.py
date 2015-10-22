#!/usr/bin/env python

import os
import matplotlib.pyplot as plt
import camb_utils
from camb_utils.cambio import cbpowerspectrum
from camb_utils import example_data
# example_data =  os.path.join(os.path.split(camb_utils.__file__)[0], 'example_data')
fname = os.path.join(example_data, 'oneh_transfer_out.dat')
ps = cbpowerspectrum(transferfile=fname,
                     Omegacdm=0.3,
                     Omegab=0.05,
                     h=0.71,
                     Omeganu=0.0,
                     As=None,
                     ns=None,
                     koverh=None)

plt.loglog(ps[:, 0], ps[:, 1])
plt.show()
