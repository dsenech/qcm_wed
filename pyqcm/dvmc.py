import pyqcm
import numpy as np

def dvmc_solver():
    # printing the model data in a file
    np.savetxt('hopping.def', pyqcm.cluster_hopping_matrix(full=True), delimiter='\t', fmt='%1.6g', comments='')
    np.savetxt('interaction.def', pyqcm.interaction_matrix(), delimiter='\t', fmt='%1.6g', comments='')

    # calling the dvmc solver
    print('Here the DVMC solver is called...')

    # reading the solution
    # here, a file qmatrix.def is read. It contains information about the Hilbert space sector, the ground state energy, and the Q matrix 
    # ....
