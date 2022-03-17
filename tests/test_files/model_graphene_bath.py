from pyqcm import *
import numpy as np

# Constructing a new cluster model
ns = 2 # number of physical sites
nb = 4 # number of bath sites
no = ns+nb # total number of sites
new_cluster_model('clus', ns, nb)

########################### BATH OPERATORS ########################### 
# Defining the bath hopping operators
new_cluster_operator('clus', 'tb1', 'one-body', [
    (1, 3, -1.0),
    (2, 4, -1.0),
    (1+no, 3+no, -1.0),
    (2+no, 4+no, -1.0)
]) # note that the last two entries pertain to the SPIN DOWN part of the operator

new_cluster_operator('clus', 'tb2', 'one-body', [
    (1, 5, -1.0),
    (2, 6, -1.0),
    (1+no, 5+no, -1.0),
    (2+no, 6+no, -1.0)
])

# Defining the 'orbital energy' of the baths
new_cluster_operator('clus', 'eb1', 'one-body', [
    (3, 3, 1.0),
    (4, 4, 1.0),
    (3+no, 3+no, 1.0),
    (4+no, 4+no, 1.0)
])

new_cluster_operator('clus', 'eb2', 'one-body', [
    (5, 5, 1.0),
    (6, 6, 1.0),
    (5+no, 5+no, 1.0),
    (6+no, 6+no, 1.0)
])

######################################################################

add_cluster('clus', [0,0,0], [[0,0,0], [1,0,0]]) # Adding the cluster in
lattice_model('Graphene_2', [[1,-1,0], [2,1,0]], [[1,-1,0], [2,1,0]]) # Tiling like an old scissorgrid elevator
set_basis([[1,0,0],[-0.5,np.sqrt(3)/2,0]]) # Classic Graphene basis (for simplicity and graphical purposes)

# Defining the interaction operator on BOTH bands
interaction_operator('U', band1=1, band2=1)
interaction_operator('U', band1=2, band2=2)

# Defining NN hopping terms
hopping_operator('t', [1,0,0], -1, band1=1, band2=2) # All hops here are from one band to another
hopping_operator('t', [0,1,0], -1, band1=1, band2=2)
hopping_operator('t', [-1,-1,0], -1, band1=1, band2=2)

