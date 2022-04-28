from pyqcm import *
from pyqcm.loop import *

import model_1D_4
set_global_parameter('accur_OP', 1e-3)
set_target_sectors(['R0:N4:S0'])
set_parameters("""
t = 1
U = 4
mu = 8
""")

##################################################################
# TEST UNITAIRE

def test_fixed_density_loop():

    banner('testing fixed_density_loop()', c='#', skip=1)

    def F():
        new_model_instance()
        averages()

    fixed_density_loop(
        5,  # starting  value of mu
        1.1, # target density
        kappa=1.0,
        maxdmu=0.5,  # maximum change in mu
        func=F,
        loop_param='U', 
        loop_values=np.arange(6, 4, -0.2),
        dens_tol=0.002,
        dir='',
        measure=None
    )


test_fixed_density_loop()

