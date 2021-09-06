from pyqcm import *
from pyqcm.loop import *

import model_2x2_C2
set_global_parameter('accur_OP', 1e-3)

##################################################################
# TEST UNITAIRE

def test_fixed_density_loop():

    banner('testing fixed_density_loop()', c='#', skip=1)

    def F():
        new_model_instance()
        averages()

    fixed_density_loop(
        5,  # starting  value of mu
        1.05, # target density
        kappa=1.0,
        maxdmu=0.2,  # maximum change in mu
        func=F,
        loop_param='U', 
        loop_values=np.arange(6, 0, -0.1),
        # var_param='M_1',
        dens_tol=0.002,
        dir='',
        measure=None
    )


test_fixed_density_loop()

