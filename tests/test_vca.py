from pyqcm import *
from pyqcm.vca import *
import matplotlib.pyplot as plt

# ax = plt.gca()
ax = None
F = None

# print(dir(pyqcm.vca)); exit()
#-----------------------------------------------------------------

# plt.show()
##################################################################
# TEST UNITAIRE

def test_vca():
    x = None
    import model_2x2_C2

    F = 'test_vca.pdf'

    banner('testing vca()', c='#', skip=1); vca(names=['M_1', 't_1'], accur=[5e-4, 5e-4], max=[10,10], NR=False)

    banner('testing vca() with Newton-Raphson', c='#', skip=1); vca(names=['M_1', 't_1'], accur=[5e-4, 5e-4], max=[10,10], NR=True)

    banner('testing vca() with explicit starting values', c='#', skip=1); vca(names=['M_1', 't_1'], start=[0.1, 1.1], accur=[5e-4, 5e-4], max=[10,10], NR=False)

    banner('testing plot_sef()', c='#', skip=1); plot_sef('M_1', np.arange(1e-9, 0.3, 0.02), accur_SEF=1e-4, show=True)

    banner('testing plot_GS_energy()', c='#', skip=1); plot_GS_energy('M_1', np.arange(1e-9, 0.3, 0.02))

    banner('testing transition_line()', c='#', skip=1)
    set_parameter('U', 1)
    set_parameter('t2', -0.5)
    transition_line('M_1', 'U', np.arange(1, 8.01, 0.25), 'mu', [0.8, 1], delta=1, verb=True)


test_vca()


# transition
# transition_line
# vca_min