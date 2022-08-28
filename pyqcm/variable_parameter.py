import numpy as np
import pyqcm

class variable_parameter:
    """This class contains the elements needed to perform computations with parameters that depend on some observable (typically the density). The values of the parameters are adjusted as a function of the observable until some self-consistency is reached.

    attributes:
        - param (str) : name of the parameter that depends on the observable
        - F : function that expresses the parameter as a function of the observable (supplied by the user)
        - obs (str) : name of the observable (a parameter). Default is 'mu'
        - lattice (boolean) : True if lattice averages are used
        - diff : difference between successive values of the parameter
        - ave : current value of the observable
        - accur : desired accuracy on the parameter

    """

    size0 = 0

    def __init__(self, param, F, obs='mu', lattice=True, accur=1e-3):
        """

        :param str param: name of the original parameter
        :param F: function that returns the value of the parameter as a function of the observable
        :param str obs: name of the parameter whose average is the observable
        :param boolean lattice: if True, the lattice average is used, otherwise the cluster average
        :param float accur: required accuracy of the self-consistent procedure

        """

        self.param = param
        self.value = 0.0
        self.F = F
        self.obs = obs
        self.lattice = lattice
        self.diff = 1e6
        self.ave = 0
        self.accur = accur
        self.iter = 0    

    def update(self, pr=False):
        """Updates the value of the parameter based on the value of the observable
        
        :param boolean pr: if True, progress is printed on the screen

        """

        if not self.lattice:
            self.ave = pyqcm.cluster_averages()[self.obs][0]
        else:
            self.ave = pyqcm.averages()[self.obs]
        
        P0 = pyqcm.parameters()[self.param]
        P = self.F(self.ave)
        pyqcm.set_parameter(self.param, P)
        self.diff = P-P0
        self.iter += 1
        if pr:
            print('<{:s}> = {:1.3g},  {:s} --> {:1.3g}, diff = {:1.3g}'.format(self.obs, self.ave, self.param, P, self.diff))


    def converged(self):
        """Tests whether the self-consistent procedure has converged

        :return boolean: True if the mean-field procedure has converged
        
        """

        if np.abs(self.diff) < self.accur:
            return True
        else:
            return False


    def __str__(self):
        return 'parameter: '+self.param+', observable: <'+self.obs+'>'


################################################################################
def variable_parameter_self_consistency(F, var_params, maxiter=10, eps_algo=0, file='hartree.tsv'):
	"""Performs the Hartree approximation

	:param F: task to perform within the loop
	:param [class variable_parameter] var_params: list of variable_parameters (or single variable_parameters)
	:param int maxiter: maximum number of iterations in the self-consistent procedure

	"""

	global first_time

	if type(var_params) is not list:
		var_params = [var_params]

	pyqcm.banner('Variable parameter self-consistency', c='*', skip=1)
	SC_converged = False
	diff_tot = 1e6

	iter = 0
	while True:
		pyqcm.new_model_instance()
		F()
		iter += 1
		pyqcm.banner('self-consistent iteration {:d}'.format(iter), '-')
		diff_tot = 0
		SC_converged = True
		for i,C in enumerate(var_params):
			C.update(pr=True)
			diff_tot += np.abs(C.diff)
			SC_converged = SC_converged and C.converged()

		print('total difference = {:g}'.format(diff_tot))

		if SC_converged:
			pyqcm.write_summary(file)
			first_time = False
			break

		if iter > maxiter :
			raise RuntimeError('Maximum number of iterations exceeded in self-consistency for variable parameters! Aborting...')

	pyqcm.banner('Variable parameter procedure has converged', c='*')
	