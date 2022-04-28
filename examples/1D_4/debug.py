import pyqcm
from pyqcm.vca import *
# pyqcm.set_global_parameter("nosym")
import model_1D_4

pyqcm.set_global_parameter("verb_ED")

pyqcm.set_target_sectors(['R0:N4'])
pyqcm.set_parameters("""
t=1
U=4
mu = 2
hx = 0.1
""")
pyqcm.new_model_instance()
print(pyqcm.cluster_info())
I = pyqcm.new_model_instance(record=True)
I.print('record.py')
# pyqcm.print_averages(pyqcm.averages(['t', 'U']))
# plot_sef('t_1', [0,0.1,0.2])
pyqcm.first_time = True
for d in [0, 0.1, 0.2]:
    pyqcm.set_parameter('hx', d)
    pyqcm.new_model_instance()
    pyqcm.averages()
    pyqcm.write_summary('ave.tsv')



