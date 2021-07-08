import pyqcm
from pyqcm.spectral import *
from pyqcm.vca import *
import model_1D_4

pyqcm.set_global_parameter("verbose", 7)
sec = 'R0:N4:S0'
pyqcm.set_target_sectors([sec])
pyqcm.set_parameters("""
t=1
U = 4
V = 1
mu = 2
""")

pyqcm.solver='dvmc'
pyqcm.new_model_instance()
