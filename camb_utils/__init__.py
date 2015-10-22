import os
from . import cambio


# tie in the example_data directory
here = __file__
basedir = os.path.split(here)[0]
example_data = os.path.join(basedir, 'example_data')
