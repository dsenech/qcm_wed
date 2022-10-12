import pyqcm
from pyqcm.vca import *
# pyqcm.set_global_parameter("nosym")
import model_1D_4

pyqcm.set_global_parameter("verb_ED")

pyqcm.set_target_sectors(['R0:S0'])
pyqcm.set_parameters("""
t=1
U=0
mu = 0
D=1e-9
""")

print(pyqcm.ground_state())

