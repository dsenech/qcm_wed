from pyqcm import *
from pyqcm.cdmft import cdmft
from pyqcm.loop import controlled_loop

import model_graphene_bath

# Imposing half-filling at 6 particles in cluster + bath sites and setting total spin to 0
set_target_sectors(['R0:N6:S0'])

# Simulation parameters
set_parameters("""
    U=4
    mu=0.5*U
    t=1
    t_1=1e-9
    tb1_1=0.5
    tb2_1=1*tb1_1
    eb1_1=1.0
    eb2_1=-1.0*eb1_1
""")

# Setting up the range of values of U over which to loop the CDMFT simulation 
U_start = 0.2 # starting at U=0.2 because the nearly non-interacting case is too long to calculate
U_stop = 10.1
U_step = 0.2

# Defining a function that will run a cdmft procedure within controlled_loop()
def run_cdmft():
    U = parameters()['U']
    cdmft(varia=["tb1_1", "eb1_1"], wc=10, grid_type='self') # setting the bath operators as the variationnal parameters

# Looping over values of U
controlled_loop(
    func=run_cdmft, 
    varia=["tb1_1", "eb1_1"], 
    loop_param="U", 
    loop_range=(U_start, U_stop, U_step)
    )