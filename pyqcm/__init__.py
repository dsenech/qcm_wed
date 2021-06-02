import numpy as np
import qcm
import re
import time

parameter_set_str = ''

################################################################################
# EXCEPTIONS

class OutOfBoundsError(Exception):
    pass
class TooManyIterationsError(Exception):
    pass

class VarParamMismatchError(Exception):
    pass

class MissingArgError(ValueError):
    pass

class MinimizationError(Exception):
    pass

class ParseError(Exception):
    pass

class WrongArgumentError(ValueError):
    pass

################################################################################
# CLASSES

class model:
    def __init__(self):
        self.record = ''

    def __repr__(self):
        return self.record

class model_instance:
    def __init__(self, model, label):
        self.record = ''
        self.model = model
        self.label = label

    def __repr__(self):
        return self.model.record + self.record

    def write(self, A):
        self.record += str(A)

    def print(self, filename='record.py'):
        global parameter_set
        with open(filename, 'w') as f:
            f.write('from pyqcm import *\nset_global_parameter("nosym")\n' + self.model.record + parameter_set_str + self.record)


the_model = model()
################################################################################
# Produces the git hash for the current version

git_hash = 'NA'
def get_git_hash():
    global git_hash
    if git_hash != 'NA':
        return
    import os
    import subprocess
    pyqcm_dir = os.path.dirname(os.path.realpath(__file__))
    current_dir = os.getcwd()
    os.chdir(pyqcm_dir)
    try:
        git_hash = str(subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']))
        git_hash = git_hash[2:-3]
        # git_hash = str(os.system('git rev-parse --short HEAD'))
    except:
        pass
    os.chdir(current_dir)
get_git_hash()

################################################################################
def new_cluster_model(name, n_sites, n_bath, generators=None, bath_irrep=False):
    """Initiates a new model (no operators yet)

    :param str name: name to be given to the model
    :param int n_sites: number of cluster sites
    :param int n_bath: number of bath sites
    :param [[int]] generators: symmetry generators (2D array of ints)
    :param [boolean] bath_irrep: True if bath orbitals belong to irreducible representations of the symmetry group
    :return: None

    """

    global the_model
    the_model.record += "new_cluster_model('"+name+"', "+str(n_sites)+', '+str(n_bath)+', generators='+str(generators)+', bath_irrep='+str(bath_irrep)+')\n'

    qcm.new_model(name, n_sites, n_bath, generators, bath_irrep)


################################################################################
def new_cluster_operator(name, op_name, op_type, elem):

    """creates a new operator from its matrix elements
    
    :param str name: name of the cluster model to which the operator will belong
    :param str op_name: name of the operator
    :param str op_type: type of operator ('one-body', 'anomalous', 'interaction', 'Hund', 'Heisenberg', 'X', 'Y', 'Z')
    :param [(int,int,float)] elem: array of matrix elements (list of tuples)
    :return: None
        
    """

    global the_model
    the_model.record += "new_cluster_operator('"+name+"', '"+op_name+"', '"+op_type+"', "+str(elem)+')\n'

    if op_type == 'anomalous':
        for x in elem:  
            if x[0] >= x[1] :
                print('anomalous matrix elements of ',op_name, ' must be such that row index < column index')
                exit()

    qcm.new_operator(name, op_name, op_type, elem)

################################################################################
def new_cluster_operator_complex(name, op_name, op_type, elem):

    """creates a new operator from its complex-valued matrix elements

    :param str name: name of the cluster model to which the operator will belong
    :param str op_name: name of the operator
    :param str op_type: type of operator ('one-body', 'anomalous', 'interaction', 'Hund', 'Heisenberg', 'X', 'Y', 'Z')
    :param [(int, int, complex)] elem: array of matrix elements (list of tuples)
    :return: None

    """

    global the_model
    the_model.record += "new_cluster_operator_complex('"+name+"', '"+op_name+"', '"+op_type+"', "+str(elem)+')\n'

    qcm.new_operator_complex(name, op_name, op_type, elem)


################################################################################
def susceptibility_poles(op_name, label=0):
    """computes the dynamic susceptibility of an operator

    :param str name: name of the operator
    :param int label: label of cluster model instance
    :returns [(float,float)]: array of 2-tuple (pole, residue)

    """

    return qcm.susceptibility_poles(op_name, label)

################################################################################
def susceptibility(op_name, freqs, label=0):
    """computes the dynamic susceptibility of an operator

        :param str op_name: name of the operator
        :param [complex] freqs: array of complex frequencies
        :para int label: label of cluster model instance
        :return: array of complex susceptibilities

    """

    return qcm.susceptibility(op_name, freqs, label)

################################################################################
def qmatrix(label=0):
    """Returns the Lehmann representation of the Green function

        :param int label: label of the cluster model instance
        :return: 2-tuple made of
            1. the array of M real eigenvalues, M being the number of poles in the representation
            2. a rectangular (L x M) matrix (real of complex), L being the dimension of the Green function

    """	

    return qcm.qmatrix(label)

################################################################################
def hybridization_Lehmann(label=0):
    """Returns the Lehmann representation of the hybridization function

    :param int label: label of the cluster model instance
    :return: 2-tuple made of

        1. the array of M real eigenvalues, M being the number of poles in the representation
        2. a rectangular (L x M) matrix (real of complex), L being the dimension of the Green function

    """	

    return qcm.hybridization(label)

################################################################################
def write_cluster_instance_to_file(filename, clus=0):
    """Writes the solved cluster model instance to a text file
    
    :param str filename: name of the file
    :param int clus: label of the cluster model instance
    
    :return: None

    """

    qcm.write_instance_to_file(filename, clus)

############################   WRAPPERS for qcm   ##############################
def add_cluster(name, pos, sites, ref=0):
    """Adds a cluster to the repeated unit

    :param str name:  name of the cluster model
    :param [int] pos: base position of cluster (array of ints)
    :param [[int]] sites: list of positions of sites (2D array of ints)
    :param [int] ref: label of a previous cluster (starts at 1) to which this one is entirely equivalent (0 = no equivalence)
    :return: None

    """

    global the_model
    the_model.record += "add_cluster('"+name+"', "+str(pos)+', '+str(sites)+', ref = '+str(ref)+')\n'

    qcm.add_cluster(name, pos, sites, ref)

############################   WRAPPERS for qcm   ##############################
def print_graph(name, sites):
    """prints a graphiz (dot) program for the cluster 

    :param str name:  name of the cluster model
    :param [[int]] sites: list of positions of sites (2D array of ints)
    :return: None

    """

    qcm.print_graph(name, sites)

################################################################################
def averages(label=0, print_to_file=True):
    """
    Computes the lattice averages of the operators present in the model

    :param int label:  label of the model instance
    :param boolean print_to_file: if True, appends a line in the file 'averages.tsv'
    :return {str,float}: a dict giving the values of the averages for each parameter

    """
    return qcm.averages(label, print_to_file)

################################################################################
def cluster_Green_function(cluster, z, spin_down=False, label=0):
    """Computes the cluster Green function

    :param int cluster: label of the cluster (0 to the number of clusters-1)
    :param complex z: frequency
    :param boolean spin_down: true is the spin down sector is to be computed (applies if mixing = 	4)
    :param int label:  label of the model instance (default 0)
    :return: a complex-valued matrix

    """

    return qcm.cluster_Green_function(cluster, z, spin_down, label)

################################################################################
def cluster_Green_function_average(cluster=0, spin_down=False):
    """Computes the cluster Green function average (integral over frequencies)

    :param int cluster: label of the cluster (0 to the number of clusters-1)
    :param boolean spin_down: true is the spin down sector is to be computed (applies if mixing = 	4)
    :return: a complex-valued matrix

    """
    return qcm.Green_function_average(spin_down, cluster)

################################################################################
def cluster_self_energy(cluster, z, spin_down=False, label=0):
    """Computes the cluster self-energy

    :param int cluster: label of the cluster (0 to the number of clusters -1)
    :param complex z: frequency
    :param boolean spin_down: true is the spin down sector is to be computed (applies if mixing = 	4)
    :param int label:  label of the model instance (default 0)
    :return: a complex-valued matrix

    """

    return qcm.cluster_self_energy(cluster, z, spin_down, label)

################################################################################
def cluster_info():
    """
    :return:A list of 3-tuples: (str, int, int): name of the cluster model, number of sites, dimension of the Green function
    """
    return qcm.cluster_info()


################################################################################
def cluster_hopping_matrix(clus=0, spin_down=False, label=0):
    """
    returns the one-body matrix of cluster no i for instance 'label'

    :param cluster: label of the cluster (0 to the number of clusters - 1)
    :param boolean spin_down: True is the spin down sector is to be computed (applies if mixing = 4)
    :param int label:  label of the model instance
    :return: a complex-valued matrix
    """
    return qcm.cluster_hopping_matrix(clus, spin_down, label)

################################################################################
def CPT_Green_function(z, k, spin_down=False, label=0):
    """
    computes the CPT Green function at a given frequency

    :param z: complex frequency
    :param k: single wavevector (ndarray(3) or array of wavevectors (ndarray(N,3))
    :param boolean spin_down: True is the spin down sector is to be computed (applies if mixing = 4)
    :param int label:  label of the model instance
    :return: a single or an array of complex-valued matrices
    """
    return qcm.CPT_Green_function(z, k, spin_down, label)

################################################################################
def CPT_Green_function_inverse(z, k, spin_down=False, label=0):
    """
    computes the inverse CPT Green function at a given frequency

    :param z: complex frequency
    :param k: array of wavevectors (ndarray(N,	:param)
    :param boolean spin_down: True is the spin down sector is to be computed (applies if mixing = 4)
    :param int label:  label of the model instance
    :return: a single or an array of complex-valued matrices
    """
    return qcm.CPT_Green_function_inverse(z, k, spin_down, label)


################################################################################
def dispersion(k, spin_down=False, label=0):
    """
    computes the dispersion relation for a single or an array of wavevectors

    :param wavevector k: single wavevector (ndarray(3) or array of wavevectors (ndarray(N,3))
    :param boolean spin_down: True is the spin down sector is to be computed (applies if mixing = 4)
    :param int label:  label of the model instance
    :return: a single (ndarray(d)) or an array (ndarray(N,d)) of real values (energies). d is the reduced GF dimension.

    """
    return qcm.dispersion(k, spin_down, label)

################################################################################
def dos(z, label=0):
    """
    computes the density of states at a given frequency.

    :param complex z: frequency
    :param int label:  label of the model instance
    :return: ndarray(d) of real values, d being the reduced GF dimension

    """
    return qcm.dos(z, label)


################################################################################
def Green_function_dimension():
    """
    returns the dimension of the CPT Green function matrix
    :return: int

    """
    return qcm.Green_function_dimension()

################################################################################
def cluster_Green_function_dimension(clus=0):
    """
    returns the dimension of the cluster Green function matrix
    :param int clus: label of the cluster
    :return: int

    """
    return qcm.Green_function_dimensionC(clus)

################################################################################
def Green_function_solve(label=0):
    """
    Usually, the Green function representation is computed only when needed, in a just-in-time fashion (i.e. in a lazy way).
    This forces the computation of the Green function representation for the current instance (i.e. non lazy). 

    :param int label:  label of the model instance
    :return: None

    """
    return qcm.Green_function_solve(label)

################################################################################
def ground_state():
    """
    :return: a list of pairs (float, str) of the ground state energy and sector string, for each cluster of the system

    """
    return qcm.ground_state()

################################################################################
def cluster_averages(label=0):
    """
    Computes the average and variance of all operators of the cluster model in the cluster ground state.

    :param int label: label of the cluster model instance
    :return: a dict str : (float, float) with the averages and variances as a function of operator name

    """
    return qcm.cluster_averages(label)

################################################################################
def Lehmann_Green_function(k, band = 1, spin_down=False, label=0):
    """
    computes the Lehmann representation of the periodized Green function for a set of wavevectors

    :param k: single wavevector (ndarray(3)) or array of wavevectors (ndarray(N,3))
    :param int band: band index (starts at 1)
    :param boolean spin_down: True is the spin down sector is to be computed (applies if mixing = 4)
    :param int label:  label of the model instance
    :return: a list of pairs {poles, residues}, each of poles and residues being itself a list.

    """
    return qcm.Lehmann_Green_function(k, band, spin_down, label)

################################################################################
def hybridization_function(clus, z, spin_down=False, label=0):
    """
    returns the hybridization function for cluster 'cluster' and instance 'label'

    :param int clus: label of the cluster (0 to the number of clusters-1)
    :param complex z: frequency
    :param boolean spin_down: True is the spin down sector is to be computed (applies if mixing = 4)
    :param int label:  label of the model instance
    :return: a complex-valued matrix

    """
    return qcm.hybridization_function(z, spin_down, clus, label)

################################################################################
def lattice_model(name, superlattice, lattice=None):
    """
    initiates the lattice model.

    :param str name: the name of the model
    :param [[int]] superlattice: array of integers of shape (d,3), d being the dimension
    :param [[int]] lattice: array of integers of shape (d,3), d being the dimension. If None, will be replaced by the unit lattice.
    :return: None

    """
    global the_model
    the_model.record += "lattice_model('"+name+"', "+str(superlattice)+', '+str(lattice)+')\n'

    
    qcm.lattice_model(name, superlattice, lattice)

################################################################################
def mixing():
    """
    returns the mixing state of the system:

    * 0 -- normal.  GF matrix is n x n, n being the number of sites
    * 1 -- anomalous. GF matrix is 2n x 2n
    * 2 -- spin-flip.  GF matrix is 2n x 2n
    * 3 -- anomalous and spin-flip (full Nambu doubling).  GF matrix is 4n x 4n
    * 4 -- up and down spins different.  GF matrix is n x n, but computed twice, with spin_down = false and true

    :return: int

    """
    return qcm.mixing()

################################################################################
def model_size():
    """
    :return: a 5-tuple:

    1. the size of the supercell
    2. the number of bands
    3. a tuple containing the sizes of each cluster
    4. a tuple containing the sizes of each cluster's bath
    5. a tuple containing the references of each cluster (label of reference cluster, from 0 to the nunber of clusters-1)

    """

    return qcm.model_size()

################################################################################
def momentum_profile(name, k, label=0):
    """
    computes the momentum-resolved average of an operator

    :param str name: name of the lattice operator
    :param k: array of wavevectors (ndarray(N,	:param)
    :param int label:  label of the model instance
    :return: an array of values

    """

    return qcm.momentum_profile(name, k, label)

################################################################################
def new_model_instance(label=0, record=False):
    """
    Creates a new instance of the lattice model, with values associated to terms of the Hamiltonian.

    :param int label:  label of the model instance
    :param boolean record: if True, keeps a str-valued record of the instance in memory, containing all the data necessary to read the instance back without solving.
    :return (class model_instance): an instance of the class `model_instance`

    """

    the_instance = model_instance(the_model, label)
    qcm.new_model_instance(label)

    if record:
        params = parameter_set(True)
        for x in params:
            if params[x][1] is None:
                the_instance.record += 'set_parameter("'+str(x)+'", '+str(params[x][0])+')\n'
            
        the_instance.record += '\nnew_model_instance('+str(label)+')\n'
        mod = model_size()
        nclus = len(mod[2])
        the_instance.record += '\nsolution=[None]*'+str(nclus)+'\n'
        for i in range(nclus):
            if mod[4][i] != i:
                continue
            clabel = label*nclus+i
            the_instance.record += '\n#--------------------- cluster no '+str(i+1)+' -----------------\n'
            the_instance.record += 'solution['+str(i)+'] = """\n' + qcm.write_instance(clabel) + '\n"""\n'

        for i in range(nclus):
            if mod[4][i] != i:
                continue
            clabel = label*nclus+i
            the_instance.record += 'read_cluster_model_instance(solution['+str(i)+'], '+str(clabel)+')\n'
            

    return the_instance


################################################################################
def new_cluster_model_instance(name, values, sec, label=0):
    """Initiates a new instance of the cluster model

    :param str name: name of the cluster model
    :param {str,float} values: values of the operators
    :param str sec: target Hilbert space sectors
    :param int label: label of model_instance

    :return: None

    """

    qcm.new_model_instanceC(name, values, sec, label)

################################################################################
def read_cluster_model_instance(S, label=0):
    """reads the solution from a string 

    :param str S: long string containing the solution 
    :param int label: label of model_instance

    :return: None

    """
    qcm.read_instance(S, label)

################################################################################
def set_parameters(params, dump=True):
    """
    Defines a new set of parameters, including dependencies

      :param tuple/str params: the values/dependence of the parameters (array of 2- or 3-tuples), or string containing syntax
      :param boolean dump: if True, sets the global str parameter_set_str tothe value
    
    :return [tuple]: list of tuples of the form (str, float) or (str, float, str). The first form gives the parameter name and its value. The second gives the parameter name, a multiplier and the name of the reference parameter. See the documentation on the hierarchy of parameters.
    
    """
    global parameter_set_str


    if type(params) is str:
        elems = []
        param_set = {}
        for p in re.split('[,;\n]', params):
            if len(p.strip()) == 0:
                continue
            if p[0] == '#':
                continue
            s = p.split('=')
            if len(s) != 2:
                raise ParseError(p)
            s2 = s[1].split('*')
            param_name = s[0].strip()
            if param_name in param_set:
                raise ParseError('parameter '+param_name+' has already been assigned!')
            if len(s2)==1:  # no dependence
                elem = (param_name, float(s[1].strip()))       
            elif len(s2)==2:
                elem = (param_name, float(s2[0].strip()), s2[1].strip())
            else:
                raise ParseError(p)
            elems.append(tuple(elem))
        if dump:
            parameter_set_str = 'set_parameters("""\n'+str(params)+'""")\n'
        qcm.set_parameters(elems)
        return elems
    else:	
        qcm.set_parameters(params)
        return params

################################################################################
def set_target_sectors(sec):
    """Define the Hilbert space sectors in which to look for the ground state

    :param [str] sec: the target sectors

    :return: None

    """
    
    global the_model
    the_model.record += """
try:
    import model_extra
except:
    pass		
"""
    the_model.record += 'set_target_sectors('+str(sec)+')\n'

    qcm.set_target_sectors(sec)

################################################################################
def parameters(label=0):
    """
    returns the values of the parameters in a given instance

    :param int label:  label of the model instance
    :return: a dict {string,float}

    """
    return qcm.parameters(label)

################################################################################
def cluster_parameters(label=0):
    """
    returns the values of the cluster parameters in a given instance, as well as the cluster model name

    :param int label:  label of the cluster model instance
    :return: a tuple:  dict{string,float}, str

    """
    return qcm.parametersC(label)

################################################################################
def parameter_set(opt='all'):
    """
    returns the content of the parameter set

    :param str opt: governs the action of the function
    :return: depends on opt

    if opt = 'all', all parameters as a dictionary {str,(float, str, float)}. The three components are 

    (1) the value of the parameter, 
    (2) the name of its overlord (or None), 
    (3) the multiplier by which its value is obtained from that of the overlord.

    if opt = 'independent', returns only the independent parameters, as a dictionary {str,float}
    if opt = 'report', returns a string with parameter values and dependencies.

    """
    P = qcm.parameter_set()
    if opt == 'independent':
        P2 = {}
        for x in P:
            if P[x][1] == None:
                P2[x] = P[x][0]
        return P2
    elif opt == 'report':
        rep = ''
        for x in P:
            if P[x][1] == None:
                rep += x + ' = ' + str(P[x][0]) + '\n'
            else:
                rep += x + ' = ' + str(P[x][2]) + ' x ' + P[x][1] + '\n'
        return rep
    else:
        return P
        

################################################################################
def periodized_Green_function(z, k, spin_down=False, label=0):
    """
    computes the periodized Green function at a given frequency and wavevectors

    :param complex z: frequency
    :param k: single wavevector (ndarray(3) or array of wavevectors (ndarray(N,3))
    :param boolean spin_down: true is the spin down sector is to be computed (applies if mixing = 4)
    :param int label:  label of the model instance
    :return: a single (d,d) or an array (N,d,d) of complex-valued matrices. d is the reduced GF dimension.

    """
    return qcm.periodized_Green_function(z, k, spin_down, label)

################################################################################
def band_Green_function(z, k, spin_down=False, label=0):
    """
    computes the periodized Green function at a given frequency and wavevectors, in the band basis (defined
    in the noninteracting model). It only differs from the periodized Green function in multi-band models.

    :param complex z: frequency
    :param k: single wavevector (ndarray(3) or array of wavevectors (ndarray(N,3))
    :param boolean spin_down: true is the spin down sector is to be computed (applies if mixing = 4)
    :param int label:  label of the model instance
    :return: a single (d,d) or an array (N,d,d) of complex-valued matrices. d is the reduced GF dimension.

    """
    return qcm.band_Green_function(z, k, spin_down, label)

################################################################################
def periodized_Green_function_element(r, c, z, k, spin_down=False, label=0):
    """
    computes the element (r,c) of the periodized Green function at a given frequency and wavevectors (starts at 0)

    :param int r: a row index (starts at 0)
    :param int c: a column index (starts at 0)
    :param complex z: frequency
    :param k: array of wavevectors (ndarray(N,3))
    :param boolean spin_down: true is the spin down sector is to be computed (applies if mixing = 4)
    :param int label:  label of the model instance
    :return: a vector of complex numbers

    """
    return qcm.periodized_Green_function_element(r, c, z, k, spin_down, label)


################################################################################
def Potthoff_functional(label=0):
    """
    computes the Potthoff functional for a given instance

    :param int label: label of the model instance
    :return: the value of the self-energy functional

    """
    return qcm.Potthoff_functional(label)

################################################################################
def properties(label=0):
    """
    Returns two strings of properties of a model instance
    
    :param int label:  label of the model instance
    :return: a pair of strings (the description line and the data line).

    """
    des, data = qcm.properties(label)
    des += 'githash_pyqcm\t'
    data += git_hash + '\t'
    return des, data

################################################################################
def print_options(opt=0):
    """Prints the list of global options and parameters on the screen

      :param int opt: 0 -> prints to screen. 1 -> prints to latex. 2 -> prints to RST

    """
    return qcm.print_options(opt)

################################################################################
def print_model(filename, **kwargs):
    """Prints a description of the model into a file

    :param str filename: name of the file

    :return: None
 
   :Keyword Arguments:

        * asy_operators (``str``) -- true if asymptote files are produced for each operator
        * asy_labels (``str``) -- true if sites labels are indicated in the asymptote program
        * asy_band (``str``) -- true if band labels are indicated in the asymptote program
        * asy_neighbors (``str``) -- true if neighbors are drawn in the asymptote program
        * asy_working_basis (``str``) -- true if the working basis is used instead of the physical basis

    """
    qcm.print_model(filename, **kwargs)

################################################################################
def projected_Green_function(z, spin_down=False, label=0):
    """
    computes the projected Green function at a given frequency, as used in CDMFT.

    :param complex z: frequency
    :param boolean spin_down: true is the spin down sector is to be computed (applies if mixing = 4)
    :param int label:  label of the model instance
    :return: the projected Green function matrix (d x d), d being the dimension of the CPT Green function.

    """
    return qcm.projected_Green_function(z, spin_down, label)

################################################################################<<<<<<<<< eliminer
def read_model(filename):
    """
    Reads the definition of the model from a file

    :param file: name of the file, the same as that of the model, i.e., without the '.model' suffix.
    :return: None

    """
    qcm.read_model(filename)

################################################################################
def reduced_Green_function_dimension():
    """
    returns the dimension of the reduced Green function, i.e. a simple multiple of the
    number of bands n, depending on the mixing state: n, 2n or 4n, and the number of bands

    """

    d = qcm.reduced_Green_function_dimension()
    nb = d
    mix = mixing()
    if mix == 1 or mix == 2:
        nb = d//2
    if mix == 3:
        nb = d//4
    return d, nb


################################################################################
def self_energy(z, k, spin_down=False, label=0):
    """
    computes the self-energy associated with the periodized Green function at a given frequency and wavevectors

    :param complex z: frequency
    :param k: single wavevector (ndarray(3) or array of wavevectors (ndarray(N,3))
    :param boolean spin_down: True is the spin down sector is to be computed (applies if mixing = 4)
    :param int label:  label of the model instance 
    :return: a single (d,d) or an array (N,d,d) of complex-valued matrices. d is the reduced GF dimension.

    """
    return qcm.self_energy(z, k, spin_down, label)

################################################################################
def set_global_parameter(name, value=None):
    """
    sets the value of a global parameter. 
    If the global parameter is Boolean, then value is True by default and should not be specified.

    :param str name: name of the global option
    :param value: value of that option (None, int, float or str)
    :return: None

    """
    if value == None:
        return qcm.set_global_parameter(name)
    else:
        return qcm.set_global_parameter(name, value)

################################################################################
def set_parameter(name, value, pr=False):
    """
    sets the value of a parameter within a parameter_set

    :param str name: name of the parameter
    :param float value: its value
    :return: None

    """
    if pr:
        print('-----> ', name, ' = ', value)
    qcm.set_parameter(name, value)

################################################################################
def spatial_dimension():
    """
    returns the spatial dimension (0, 1, 2, 3) of the model

    """
    return qcm.spatial_dimension()

################################################################################
def spectral_average(name, z, label=0):
    """
    returns the contribution of a frequency to the average of an operator

    :param name: name of the operator
    :param z: complex frequency
    :param int label:  label of the model instance
    :return: float

    """
    return qcm.spectral_average(name, z, label)

################################################################################
def variational_parameters():
    """
     
    :return: a list of the names of the variational parameters

    """
    return qcm.variational_parameters()

################################################################################
def set_basis(B):
    """
    
    :param B: the basis (a (D x 3) real matrix)
    :return: None

    """

    global the_model
    the_model.record += "set_basis("+str(B)+')\n'

    qcm.set_basis(B)

################################################################################
def interaction_operator(name, **kwargs):
    """
    Defines an interaction operator of type Hubbard, Hund, Heisenberg or X, Y, Z

    :param str name: name of the operator

   :Keyword Arguments:

        * link (*[int]*): 3-component integer vector, (0,0,0) by default
        * amplitude (*float*): amplitude multiplier
        * band1 (*int*): Band label of first index (1 by default)
        * band2 (*int*): Band label of second index (1 by default)
        * type (*str*): one of 'Hubbard', 'Heisenberg', 'Hund', 'X', 'Y', 'Z'

    :return: None

    """

    global the_model
    the_model.record += "interaction_operator('"+name+"'"
    for x in kwargs:
        if type(kwargs[x]) is str:
            the_model.record += ', '+x+"='"+kwargs[x]+"'"
        else:	
            the_model.record += ', '+x+'='+str(kwargs[x])
    the_model.record += ')\n'	

    qcm.interaction_operator(name, **kwargs)


################################################################################
def hopping_operator(name, link, amplitude, **kwargs):
    """Defines a hopping term or, more generally, a one-body operator

    :param str name: name of operator
    :param [int] link: bond vector (3-component integer array)
    :param float amplitude: hopping amplitude multiplier
    
    :Keyword Arguments:

        * band1 (*int*) -- Band label of first index (1 by default)
        * band2 (*int*) -- Band label of second index (1 by default)
        * tau (*int*) -- specifies the tau Pauli matrix  (0,1,2,3)
        * sigma (*int*) -- specifies the sigma Pauli matrix  (0,1,2,3)
  
    :return: None

    """

    global the_model
    the_model.record += "hopping_operator('"+name+"', "+str(link)+', '+str(amplitude)
    for x in kwargs:
        if type(kwargs[x]) is str:
            the_model.record += ', '+x+"='"+kwargs[x]+"'"
        else:	
            the_model.record += ', '+x+'='+str(kwargs[x])
    the_model.record += ')\n'	

    qcm.hopping_operator(name, link, amplitude, **kwargs)

################################################################################
def anomalous_operator(name, link, amplitude, **kwargs):
    """Defines an anomalous operator

      :param str name: name of operator
      :param [int] link: bond vector (3-component integer array)
      :param complex amplitude: pairing multiplier

    :Keyword Arguments:

        * band1 (*int*) -- Band label of first index (1 by default)
        * band2 (*int*) -- Band label of second index (1 by default)
        * type (*str*) -- one of 'singlet' (default), 'dz', 'dy', 'dx'
  
    :return: None

    """

    global the_model
    the_model.record += "anomalous_operator('"+name+"', "+str(link)+', '+str(amplitude)
    for x in kwargs:
        if type(kwargs[x]) is str:
            the_model.record += ', '+x+"='"+kwargs[x]+"'"
        else:	
            the_model.record += ', '+x+'='+str(kwargs[x])
    the_model.record += ')\n'	

    qcm.anomalous_operator(name, link, amplitude, **kwargs)

################################################################################
def explicit_operator(name, elem, **kwargs):
    """
    Defines an explicit operator

    :param str name: name of operator
    :param [(int,int,complex)] elem: matrix elements

    :Keyword Arguments:

        * tau (*int*) -- specifies the tau Pauli matrix  (0,1,2,3)
        * sigma (*int*) -- specifies the sigma Pauli matrix  (0,1,2,3)
        * type (*str*) -- one of 'one-body' [default], 'singlet', 'dz', 'dy', 'dx', 'Hubbard', 'Hund', 'Heisenberg', 'X', 'Y', 'Z'

    :return: None

    """
    global the_model
    the_model.record += "explicit_operator('"+name+"', "+str(elem)
    for x in kwargs:
        if type(kwargs[x]) is str:
            the_model.record += ', '+x+"='"+kwargs[x]+"'"
        else:	
            the_model.record += ', '+x+'='+str(kwargs[x])
    the_model.record += ')\n'	

    qcm.explicit_operator(name, elem, **kwargs)

################################################################################
def density_wave(name, t, Q, **kwargs):
    """
    Defines a density wave

    :param str name: name of operator
    :param str t: type of density-wave -- one of 'Z', 'X', 'Y', 'cdw', 'singlet', 'dz', 'dy', 'dx'
    :param wavevector Q:  wavevector of the density wave (in multiple of :math:`pi`)

    :Keyword Arguments:

        * link (*[int]*) -- bond vector, for bond density waves
        * amplitude (*complex*) -- amplitude multiplier
        * band (*int*) -- Band label (1 by default)
        * phase (*float*) -- real phase (as a multiple of :math:`pi`)

    :return: None

    """

    global the_model
    the_model.record += "density_wave('"+name+"', '"+str(t)+"', "+str(Q)
    for x in kwargs:
        if type(kwargs[x]) is str:
            the_model.record += ', '+x+"='"+kwargs[x]+"'"
        else:	
            the_model.record += ', '+x+'='+str(kwargs[x])
    the_model.record += ')\n'	

    qcm.density_wave(name, t, Q, **kwargs)

################################################################################
def site_and_bond_profile():
    """
    Computes the site and bond profiles in all clusters of the repeated unit

    :return: A pair of ndarrays

    site profile -- the components are 
    x y z n Sx Sy Sz psi.real psi.imag

    bond profile -- the components are  
    x1 y1 z1 x2 y2 z2 b0 bx by bz d0.real dx.real dy.real dz.real d0.imag dx.imag dy.imag dz.imag

    """
    return qcm.site_and_bond_profile()

################################################################################
def V_matrix(z, k, spin_down=False, label=0):
    """
    Computes the matrix :math:`V=G_0^{-1}-G^{c-1}_0` at a given frequency and wavevectors, where :math:`G_0` is the noninteracting Green function on the infinite lattice and :math:`G^c_0` is the noninteracting Green function on the cluster.

    :param complex z: frequency
    :param wavevector k: wavevector (ndarray(3))
    :param boolean spin_down: True is the spin down sector is to be computed (applies if mixing = 4)
    :param int label:  label of the model instance
    :return: a single (d,d) or an array (N,d,d) of complex-valued matrices. d is the reduced GF dimension.

    """
    return qcm.V_matrix(z, k, spin_down, label)

################################################################################
def tk(k, spin_down=False, label=0):
    """
    computes the k-dependent one-body matrix of the lattice model

    :param k: single wavevector (ndarray(3) or array of wavevectors (ndarray(N,3))
    :param boolean spin_down: True is the spin down sector is to be computed (applies if mixing = 4)
    :param int label:  label of the model instance
    :return: a single or an array of complex-valued matrices
    """
    return qcm.tk(k, spin_down, label)

################################################################################
def QP_weight(k, eta=0.01, band=1, spin_down=False, label=0):
    """
    computes the k-dependent quasi-particle weight from the self-energy derived from the periodized Green function

    :param k: single wavevector (ndarray(3) or array of wavevectors (ndarray(N,3))
    :param float eta: increment in the imaginary axis direction used to computed the derivative of the self-energy
    :param int band: band index (starts at 1)
    :param boolean spin_down: True is the spin down sector is to be computed (applies if mixing = 4)
    :param int label:  label of the model instance
    :return: a single float or an array of floats, depending on the shape of k

    """
    sigma1 = qcm.self_energy(-eta*1j, k, spin_down, label)
    sigma2 = qcm.self_energy(eta*1j, k, spin_down, label)
    if len(sigma1.shape) == 3:
        Z = (sigma1[:,band-1,band-1].imag - sigma2[:,band-1,band-1].imag)/(2*eta) + np.ones(len(k))
    else:
        Z = (sigma1[band-1,band-1].imag - sigma2[band-1,band-1].imag)/(2*eta) + 1.0
    Z = 1.0/Z
    return Z

################################################################################
def cluster_QP_weight(cluster=0, eta=0.01, band=1, spin_down=False, label=0):
    """
    computes the cluster quasi-particle weight from the cluster self-energy

    :param int cluster: cluster label (starts at 0)
    :param float eta: increment in the imaginary axis direction used to computed the derivative of the self-energy
    :param int band: band index (starts at 1)
    :param boolean spin_down: True is the spin down sector is to be computed (applies if mixing = 4)
    :param int label:  label of the model instance
    :return: a float

    """
    sigma1 = cluster_self_energy(cluster, -eta*1j, spin_down, label)
    sigma2 = cluster_self_energy(cluster, eta*1j, spin_down, label)
    Z = (sigma1[band-1,band-1].imag - sigma2[band-1,band-1].imag)/(2*eta) + 1.0
    Z = 1.0/Z
    return Z
    

################################################################################
def spin_spectral_function(freq, k, band=1, label=0):
    """
    computes the k-dependent spin-resolved spectral function

    :param freq: complex freqency
    :param k: single wavevector (ndarray(3) or array of wavevectors (ndarray(N,3))
    :param int band: band index (starts at 1)
    :param int label:  label of the model instance
    :return: depending on the shape of k, a nd.array(3) of nd.array(N,3)

    """
    mix = mixing()
    ds, nbands = reduced_Green_function_dimension()
    if mix != 2 and mix != 3:
        print('The function "spin_spectral_function()" makes sense only if spin-flip terms are present')
        exit()
    if mix == 2:
        ds //= 2
    elif mix == 3:
        ds //= 4
    
    G = periodized_Green_function(freq, k, False, label)

    if len(G.shape) == 3:
        nk = G.shape[0]
        S = np.zeros((nk,4))
        if band is None:
            for l in range(ds):
                s1 = G[:, l, l+ds]
                s2 = G[:, l+ds, l]
                S[:,1] += -(s1+s2).imag
                S[:,2] +=  (s1-s2).real
                S[:,3] += -(G[:, l, l] - G[:, l+ds, l+ds]).imag
                S[:,0] += -(G[:, l, l] + G[:, l+ds, l+ds]).imag
        else:
            l = band-1
            s1 = G[:, l, l+ds]
            s2 = G[:, l+ds, l]
            S[:,1] = -(s1+s2).imag
            S[:,2] =  (s1-s2).real
            S[:,3] = -(G[:, l, l] - G[:, l+ds, l+ds]).imag
            S[:,0] = -(G[:, l, l] + G[:, l+ds, l+ds]).imag
    else:
        S = np.zeros(4)
        if band is None:
            for l in range(ds):
                s1 = G[l, l+ds]
                s2 = G[l+ds, l]
                S[1] += -(s1+s2).imag
                S[2] +=  (s1-s2).real
                S[3] += -(G[l, l] - G[l+ds, l+ds]).imag
                S[0] += -(G[l, l] + G[l+ds, l+ds]).imag
        else:
            l = band-1
            s1 = G[l, l+ds]
            s2 = G[l+ds, l]
            S[1] = -(s1+s2).imag
            S[2] =  (s1-s2).real
            S[3] = -(G[l, l] - G[l+ds, l+ds]).imag
            S[0] = -(G[l, l] + G[l+ds, l+ds]).imag

    return 0.5*S

################################################################################
def update_bath(label=0):
    """
    updates the model parameters without creating a new instance or reseting the instance specified

    :param int label:  label of the model instance
    :return: None

    """
    qcm.update_bath(label)

################################################################################
def print_wavefunction(label=0):
    """
    prints the ground state wavefunction(s) on the screen

    :param int label:  label of the model instance
    :return: None

    """
    return qcm.print_wavefunction(label)

################################################################################
def matrix_elements(model, op):
    """
    returns the type and matrix elements defining a Hermitian operator

    :param str model: name of the cluster model
    :param str op: name of the operator
    :return: a tuple (typ, elem)

    """
    return qcm.matrix_elements(model, op)



################################################################################


################################################################################
# GENERIC FUNCTIONS
################################################################################
def banner(s, c='-', skip=0):
    if skip:
        print('\n'*(skip-1))
    n = len(s)
    if n > 72 :
        print(c[0]*80, flush=True)
        print(s, flush=True)
        print(c[0]*80, flush=True)
    else:
        m = (80 - n - 2) // 2
        print(c[0]*m, s, c[0]*m, flush=True)
    if skip:
        print('\n'*(skip-1))



def print_averages(ave):
    """
    Prints the averages nicely on the screen

    :param dict ave: the dictionary produced by qcm.averages()

    """
    for x in ave:
        print('   <{:s}> = {:f}'.format(x, ave[x]))

def print_cluster_averages(ave):
    """
    Prints the averages nicely on the screen

    :param dict ave: the dictionary produced by qcm.averages()

    """
    for x in ave:
        print('   <{:s}> = {:f}'.format(x, ave[x][0]))
        

def print_parameters(P):
    """
    Prints the parameters nicely on the screen
        
    :param dict P: the dictionary of parameters
    
    """
    for x in P:
        print(x, ' = ', P[x])


def write_summary(f, first=False, suppl_descr=None, suppl_values=None):
    fout = open(f, 'a')
    des, val = properties()
    if suppl_values != None:
        val += suppl_values
    val += time.strftime("%Y-%m-%d@%H:%M", time.localtime())
    if first:
        if suppl_descr != None:
            des += suppl_descr
        des += 'time'    
        fout.write('\n\n' + des + '\n')
    fout.write(val + '\n')
    fout.close()


################################################################################
def __wavevector_line(k1, k2, n=32):
    """
    Builds a wavevector path and associated tick marks for a straight path between k1 and k2 with n points
    
    :param (float) k1 : starting wavevector
    :param (float) k2 : ending wavevector
    :param int n: number of wavevectors per segment
    :returns tuple: 1) a ndarray of wavevectors 2) a list of tick positions 3) a list of tick strings

    """
    ticks = np.array([0, n + 1])
    tick_labels = [str(k1), str(k2)]
    k1 = np.array(k1)
    k2 = np.array(k2)
    k = np.zeros((n + 1, 3))
    n1 = 1.0/n
    for i in range(n + 1):
        k[i,:] = n1*((n-i)*k1 + i*k2)
    return 0.5 * k, ticks, tick_labels
    

################################################################################
def wavevector_path(n=32, shape='triangle'):
    """
    Builds a wavevector path and associated tick marks
    
    :param int n: number of wavevectors per segment
    :param str shape: the geometry of the path, one of: line, halfline, triangle, graphene, graphene2, diagonal, cubic, cubic2, tetragonal, tetragonal2  OR a tuple with two wavevectors for a straight path between the two
    :returns tuple: 1) a ndarray of wavevectors 2) a list of tick positions 3) a list of tick strings

    """
    if type(shape) is tuple:
        return __wavevector_line(shape[0], shape[1], n=32)

    elif (shape == 'triangle'):
        k = np.zeros((3 * n + 1, 3))
        for i in range(n):
            k[i, 0] = i / n
        for i in range(n):
            k[i + n, 0] = 1.0
            k[i + n, 1] = i / n
        for i in range(n):
            k[i + 2 * n, 0] = 1.0 - i / n
            k[i + 2 * n, 1] = 1.0 - i / n
        k[-1] = k[0]
        ticks = np.array([0, n, 2 * n, 3 * n + 1])
        tick_labels = [r'$(0,0)$', r'$(\pi,0)$', r'$(\pi,\pi)$', r'$(0,0)$']
    elif shape == 'line':
        k = np.zeros((2 * n + 1, 3))
        for i in range(2 * n + 1):
            k[i, 0] = (i - n) / n
        ticks = np.array([0, n, 2 * n + 1])
        tick_labels = [r'$-\pi$', r'$0$', r'$\pi$']
    elif shape == 'diagonal':
        k = np.zeros((n + 1, 3))
        for i in range(n + 1):
            k[i, 0] = i / n
            k[i, 1] = i / n
        ticks = np.array([0, n / 2, n + 1])
        tick_labels = [r'$\Gamma$', r'$(\pi/2,\pi/2)$', r'$M$']
    elif shape == 'halfline':
        k = np.zeros((n + 1, 3))
        for i in range(n + 1):
            k[i, 0] = i / n
        ticks = np.array([0, n / 2, n + 1])
        tick_labels = [r'$0$', r'$\pi/2$', r'$\pi$']
    elif shape == 'graphene':  # honeycomb lattice (2D) gamma-M-K'-gamma
        k = np.zeros((5 * n // 2 + 1, 3))
        for i in range(n):
            k[i, 0] = i * 0.6666667 / n
        for i in range(1 + n // 2):
            k[n + i, 0] = 0.66666667
            k[n + i, 1] = -4 * i * 0.19245 / n
        for i in range(n):
            k[-i - 1, 0] = i * 0.6666667 / n
            k[-i - 1, 1] = -2 * i * 0.19245 / n
        ticks = np.array([0, n, 3 * n // 2, 5 * n // 2])
        tick_labels = [r'$\Gamma$', r'$M$', r'$K^\prime$', r'$\Gamma$']
    elif shape == 'graphene2':  # honeycomb lattice (2D) M-gamma-K-K'
        k = np.zeros((3*n + 1, 3))
        for i in range(n):
            k[i, 0] = (n-i) * 0.6666667 / n
            k[i, 1] = -2 * (n-i) * 0.19245 / n
        for i in range(n):
            k[i+n, 0] = i * 0.6666667 / n
        for i in range(n+1):
            k[2*n+i, 0] = 0.66666667
            k[2*n+i, 1] = -2 * i * 0.19245 / n
        ticks = np.array([0, n, 2*n, 3*n])
        tick_labels = [r'$K$', r'$\Gamma$', r'$M$', r'$K^\prime$']
    elif shape == 'tri':  # triangular lattice (2D) gamma-M-K'-gamma
        k = np.zeros((5 * n // 2 + 1, 3))
        sq3 = np.sqrt(3.0)
        for i in range(n):
            k[i,1] = i * 0.6666667*sq3 / n
        for i in range(1 + n // 2):
            k[n + i, 1] = 0.66666667*sq3
            k[n + i, 0] = -4 * i * 0.19245*sq3 / n
        for i in range(n):
            k[-i - 1, 1] = i * 0.6666667*sq3 / n
            k[-i - 1, 0] = -2 * i * 0.19245*sq3 / n
        ticks = np.array([0, n, 3 * n // 2, 5 * n // 2])
        tick_labels = [r'$\Gamma$', r'$M$', r'$K$', r'$\Gamma$']
    elif shape == 'cubic':  # cubic lattice (100)-(000)-(111)-(110)-(000)
        k = np.zeros((4*n+1, 3))
        for i in range(n):
            k[i, 0] = 1.0 - i / n
        for i in range(n):
            k[i + n, 0] = i / n
            k[i + n, 1] = i / n
            k[i + n, 2] = i / n
        for i in range(n):
            k[i + 2*n, 0] = 1
            k[i + 2*n, 1] = 1
            k[i + 2*n, 2] = 1.0 - i / n
        for i in range(n):
            k[i + 3*n, 0] = 1.0 - i / n
            k[i + 3*n, 1] = 1.0 - i / n
        ticks = np.array([0, n, 2 * n, 3 * n, 4 * n + 1])
        tick_labels = [r'$(\pi,0,0)$', r'$(0,0,0)$', r'$(\pi,\pi,\pi)$', r'$(\pi,\pi,0)$', r'$(0,0,0)$']
    elif shape == 'tetragonal':  # tetragonal lattice (000)-(200)-(101/2)-(111/2)-(110)
        k = np.zeros((4*n+1, 3))
        for i in range(n):
            k[i, 0] = 2 * i / n
            k[i, 1] = 0
            k[i, 2] = 0
        for i in range(n):
            k[i + n, 0] = 2 - i / (n)
            k[i + n, 1] = 0
            k[i + n, 2] = i / (2 * n)
        for i in range(n):
            k[i + 2*n, 0] = 1 
            k[i + 2*n, 1] = i / n
            k[i + 2*n, 2] = 1 / 2
        for i in range(n):
            k[i + 3*n, 0] = 1 
            k[i + 3*n, 1] = 1 
            k[i + 3*n, 2] = 1 / 2 - i / (2 * n)
        k[-1,0] = 1
        k[-1,1] = 1
        k[-1,2] = 0 
        ticks = np.array([0, n, 2 * n, 3 * n, 4 * n + 1])
        tick_labels = [r'$(0,0,0)$', r'$(2\pi,0,0)$', r'$(\pi,0,\pi/2)$', r'$(\pi,\pi,\pi/2)$', r'$(\pi,\pi,0)$']
    elif shape == 'tetragonal2':  # tetragonal lattice (000)-(100)-(1/201/4)-(1/21/21/4)-(1/21/20)
        k = np.zeros((4*n+1, 3))
        for i in range(n):
            k[i, 0] = i / n
            k[i, 1] = 0
            k[i, 2] = 0
        for i in range(n):
            k[i + n, 0] = 1 - i / (2 * n)
            k[i + n, 1] = 0
            k[i + n, 2] = i / (4 * n)
        for i in range(n):
            k[i + 2*n, 0] = 1 / 2
            k[i + 2*n, 1] = i / (2 * n)
            k[i + 2*n, 2] = 1 / 4
        for i in range(n):
            k[i + 3*n, 0] = 1 / 2
            k[i + 3*n, 1] = 1 / 2
            k[i + 3*n, 2] = 1 / 4 - i / (4 * n)
        ticks = np.array([0, n, 2 * n, 3 * n, 4 * n + 1])
        tick_labels = [r'$(0,0,0)$', r'$(\pi,0,0)$', r'$(\pi/2,0,\pi/4)$', r'$(\pi/2,\pi/2,\pi/4)$', r'$(\pi/2,\pi/2,0)$']
    elif shape == 'cubic2':  # cubic lattice (000)-(010)-(110)-(000)-(111)-(010)-(000)
        k = np.zeros((6*n+1, 3))
        for i in range(n):
            k[i, 1] = i / n
        for i in range(n):
            k[i + n, 0] = i / n
            k[i + n, 1] = 1.0
            k[i + n, 2] = 0
        for i in range(n):
            k[i + 2*n, 0] = 1.0 - i / n
            k[i + 2*n, 1] = 1.0 - i / n
            k[i + 2*n, 2] = 0
        for i in range(n):
            k[i + 3*n, 0] = i / n
            k[i + 3*n, 1] = i / n
            k[i + 3*n, 2] = i / n
        for i in range(n):
            k[i + 4*n, 0] = 1.0 - i / n
            k[i + 4*n, 1] = 1.0
            k[i + 4*n, 2] = 1.0 - i / n
        for i in range(n):
            k[i + 5*n, 0] = 0
            k[i + 5*n, 1] = 1.0 - i / n
            k[i + 5*n, 2] = 0
        k[-1,0] = 0
        k[-1,1] = 0
        k[-1,2] = 0 
        ticks = np.array([0, n, 2 * n, 3 * n, 4 * n, 5 * n, 6 * n + 1])
        tick_labels = [r'$(0,0,0)$', r'$(0,\pi,0)$', r'$(\pi,\pi,0)$', r'$(0,0,0)$', r'$(\pi,\pi,\pi)$' , r'$(0,\pi,0)$', r'$(0,0,0)$']
    else:
        print('-------> shape ', shape, ' unknown')
    return 0.5 * k, ticks, tick_labels


################################################################################
# produces a set of wavevectors for an mdc
def wavevector_grid(n=100, orig=[-1.0, -1.0], side=2, k_perp = 0, plane='z'):
    """Produces a set of wavevectors for a MDC

    :param int n: number of wavevectors on the side
    :param [float] orig: origin (in multiples of pi)
    :param float side: length of the side (in multiples of pi)
    :param float k_perp: momentum component in the third direction (in multiples of pi)
    :param str plane: momentum plane, 'xy'='z', 'yz'='x'='zy' or 'xz'='zx'='y'
    :return: ndarray of wavevectors (n*n x 3)
    
    """
    if spatial_dimension() < 2:
        print('"calling wavevector_grid()" makes no sense for a spatial dimension < 2 "')
        exit()

    c = np.array([0,1,2])
    if plane in ['y', 'xz', 'zx']:
        c[0] = 2
        c[1] = 0
        c[2] = 1
    elif plane in ['x', 'yz', 'zy']:
        c[0] = 1
        c[1] = 2
        c[2] = 0
    orig0 = 0.5*orig[0]
    orig1 = 0.5*orig[1]
    sidep = 0.5*side
    k = np.zeros((n * n, 3))
    step = 1.0 * sidep / (n-1)
    for i in range(n):
        for j in range(n):
            k[i + n * j, c[0]] = orig0 + i * step
            k[i + n * j, c[1]] = orig1 + j * step
            k[i + n * j, c[2]] = 0.5*k_perp
    return k


################################################################################
def read_from_file(out_file, n=0):
    """
    reads an output file for parameters

    :param str out_file: name of output file from which parameters are read
    :param int n: line number of data in output file (excluding titles)
    :return: string to be added to an enventual input file

    """

    print('reading data file ' + out_file + ', line ' + str(n + 1))
    fin = open(out_file, 'r')
    K = fin.readline().strip()
    K = K.split()
    V = ""
    for i in range(n + 1):
        V = fin.readline()
    V = V.split()
    print(str(len(V)) + ' columns read')

    out_keys = ['model', 'omega', 'githash', 'githash_ED', 'iterations', 'dist_function', 'diff_hybrid', 'temp', 'E_kin', 'E_pot', 'solver', 'time', 'mean_field_niter', 'distance']

    PR = ''
    for i in range(len(K)):
        if K[i] in out_keys:
            continue
        elif K[i][0:4] == 'ave_':
            continue
        elif K[i][0:4] == 'var_':
            continue
        elif K[i][0:3] == 'E0_':
            continue
        elif K[i][0:2] == 'N_':
            continue
        elif K[i][0:7] == 'sector_':
            continue
        else:
            PR += K[i] + ' = ' + str(V[i]) + '\n'

    return PR

def to_input_file(out_file, n=0):
    return read_from_file(out_file, n)

################################################################################
def params_from_file(out_file, n=0):
    """
    reads an output file for parameters

    :param str out_file: name of output file from which parameters are read
    :param int n: line number of data in output file (excluding titles).
    :return: a dict of (parameter, value)

    """

    print('reading data file ' + out_file + ', line ' + str(n + 1))
    fin = open(out_file, 'r')
    K = fin.readline().strip()
    K = K.split()
    V = ""
    for i in range(n + 1):
        V = fin.readline()
    V = V.split()
    print(str(len(V)) + ' columns read')

    out_keys = ['model', 'omega', 'githash', 'githash_ED', 'iterations', 'dist_function', 'diff_hybrid', 'temp', 'E_kin', 'E_pot', 'solver', 'time', 'mean_field_niter', 'distance']

    D = {}
    for i in range(len(K)):
        if K[i] in out_keys:
            continue
        elif K[i][0:4] == 'ave_':
            continue
        elif K[i][0:4] == 'var_':
            continue
        elif K[i][0:3] == 'E0_':
            continue
        elif K[i][0:2] == 'N_':
            continue
        elif K[i][0:7] == 'sector_':
            continue
        else:
            D[K[i]] = float(V[i])
    return D

################################################################################
def set_params_from_file(out_file, n=0):
    """
    reads an output file for parameters

    :param str out_file: name of output file from which parameters are read
    :param int n: line number of data in output file (excluding titles)
    :return: nothing

    """

    par = qcm.parameter_set()
    D = params_from_file(out_file, n)
    for x in par:
        if par[x][1] != None:
            continue
        if x in D:
            set_parameter(x,D[x],pr=True)

################################################################################
def parameter_string(lattice=True, CR=False):
    """
    Returns a string with the model parameters. 
    :param boolean lattice : if True, only indicates the independent lattice parameters.
    :param boolean CR : if True, put each parameter on a line.
    """
    par = qcm.parameter_set()
    S = ''
    sep = ', '
    if CR:
        sep = '\n'
    first = True
    for x in par:
        if par[x][1] != None:
            continue
        if '_' in x and lattice:
            continue
        if first is False:
            S += sep        
        S += x + '={:g}'.format(par[x][0])
        if first:
            first = False
    return S            


################################################################################
def read_from_file_legacy(filename):
    """
    reads model parameters from a text file, for legacy results
    :param str filename: name of the input text file
    """
    with open(filename, 'r') as f:
        F = f.read()

    try:
        i1 = F.find('parameters')
    except:
        print('"parameters" not found in "'+filename+' exiting')
        exit()

    data = F[i1:].splitlines()
    data = data[1:]

    elems = []
    for s in data:
        if len(s)==0: break
        if s[0]=='#': continue
        p = re.split('[ ,\t]+', s)
        if len(p) == 2:
            elem = (p[0].strip(), float(p[1].strip()))
        elif len(p) == 3:
            elem = (p[0].strip(), float(p[1].strip()), p[2].strip())
        else:
            break	
        elems.append(tuple(elem))
    print('parameters: ', elems)
    set_parameters(elems)

    try:
        i1 = F.find('target_sectors')
    except:
        print('"target_sectors" not found in "'+filename+' exiting')
        exit()

    i2 = F[i1:].find('\n\n')	
    data = F[i1:i2].splitlines()
    data = data[1:]
    elems = []
    for s in data:
        if s[0]=='#': continue
        terms = re.split('[, \t]', s)
        sec = ''
        for i in range(1, len(terms)):
            sec += terms[i].strip() + '/'
        elems.append(sec[:-1])
    print('target sectors: ', elems)
    set_target_sectors(elems)



################################################################################
def __varia_table(var, val, prefix = ''):
    s = prefix
    for i,p in enumerate(var):
        s += '{:<9} = {: .4f}\t'.format(p,val[i])
        if (i+1)%5 == 0:
            s += '\n'
            s += prefix
    return s
