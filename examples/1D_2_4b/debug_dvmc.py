from pyqcm import *
from pyqcm.cdmft import *

import model_1D_2_4b

np.set_printoptions(precision=4, linewidth=512, suppress=True)

sec = 'R0:N6:S0'
set_global_parameter('verbose',0)
set_target_sectors([sec])

set_parameters("""
    U=3.9
    mu=1
    t=1
    tb1_1=0.5
    tb2_1=0.5
    eb1_1=1
    eb2_1=-1
""")

pyqcm.solver='dvmc'
pyqcm.new_model_instance()

from pyqcm.spectral import *
spectral_function(path="line")
