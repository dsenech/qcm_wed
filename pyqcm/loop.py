import numpy as np
import pyqcm
import pyqcm.cdmft
import os

first_time=True

################################################################################
def loop_from_file(func, file):
	"""Performs a task over a set of instances defined in a file

	:param func: function (task) to perform
	:param str file: name of the data file specifying the solutions, in the same tab-separated format usually used to write solutions

	"""
	param = pyqcm.parameter_set()
	data = np.genfromtxt(file, names=True)
	print('number of data sets: ', len(data))

	independent_params = []
	for key, val in param.items():
		if val[1] is None and key in data.dtype.names:
			independent_params.append(key)

	print('independent parameters: ', independent_params)

	N = len(data)
	for i in range(N):
		print('*'*30, ' set ', i+1, ' of ', N, '*'*30)
		for x in independent_params:
			print(x, ' --> ', data[x][i])
			pyqcm.set_parameter(x, data[x][i])
		func()


################################################################################

def __predict_bad(x, y, p):
	for i in range(len(x)):
		if abs(x[i] - y[i]) / (abs(x[i]) + 0.01) > p:
			return True
	return False


def __predict_good(x, y, p):
	for i in range(len(x)):
		if abs(x[i] - y[i]) / abs((x[i]) + 0.01) > p:
			return False
	return True

################################################################################
# performs a loop for VCA or CDMFT with an external loop parameter and an optional control on an order parameter

def controlled_loop(
	func=None,
	varia=None, 
	loop_param=None, 
	loop_range=None, 
	control_func=None,
	adjust=False,
	predict=True
):
	"""Performs a controlled loop for VCA or CDMFT with a predictor
	
	:param func: a function called at each step of the loop
	:param [str] varia: names of the variational parameters
	:param str loop_param: name of the parameter looped over
	:param (float, float, float) loop_range: range of the loop (min, max, step)
	:param control_func: (optional) name of the function that controls the loop (returns boolean)
	:param boolean adjust: if True, adjusts the steps depending on control or status of convergence
	:param boolean predict: if True, uses a linear or quadratic predictor

	"""

	pyqcm.banner('controlled loop over '+loop_param, '%', skip=1)
	if loop_param is None:
		raise pyqcm.MissingArgError('The name of the control parameter (loop_param) is missing')
	if loop_range is None:
		raise pyqcm.MissingArgError('A loop range (min, max, step) must be provided')
	if (loop_range[1]-loop_range[0])*loop_range[2] < 0:
		raise ValueError('the loop range is incoherent: step has wrong sign')
	if varia is None:
		raise pyqcm.MissingArgError('variational parameters (varia) must be specified')
	if type(varia) != list:
		raise TypeError('the argument "varia" of controlled_loop must be a list')

	nvar = len(varia)  # number of variational parameters
	sol = np.empty((4, nvar))  # current solution + 3 previous solutions
	prm = np.empty(4)  # current loop parameter + 3 previous values

	prm[0] = loop_range[0] - loop_range[2]
	step = loop_range[2]
	first = True
	loop_counter = 0
	retry = False
	nretry = 0
	max_retry = 4
	start = [0]*nvar
	while True:  # main loop here
		if retry:
			step = step / 2
			retry = False
			nretry += 1
			if nretry > max_retry:
				break
			print('retrying with half the step')
		else:
			prm = np.roll(prm, 1)  # cycles the values of the loop parameter
			sol = np.roll(sol, 1, 0)  # cycles the solutions

		prm[0] = prm[1] + step  # updates the loop parameter
		pyqcm.set_parameter(loop_param, prm[0])
		if loop_range[1] >= loop_range[0] and prm[0] > loop_range[1]:
			break
		if loop_range[1] <= loop_range[0] and prm[0] < loop_range[1]:
			break
	
		# predicting the starting value from previous solutions (do nothing if loop_counter=0)
		fit_type = ''
		if loop_counter == 1 or not predict:
			start = sol[1, :]  # the starting point is the previous value in the loop
		elif loop_counter > 1:
			if loop_counter == 2:
				fit_type = ' (linear fit)'
				x = prm[1:3]  # the two previous values of the loop parameter
				y = sol[1:3, :]  # the two previous values of the loop parameter
				pol = np.polyfit(x, y, 1)
			else:
				fit_type = ' (quadratic fit)'
				x = prm[1:4]  # the three previous values of the loop parameter
				y = sol[1:4, :]  # the three previous values of the loop parameter
				pol = np.polyfit(x, y, 2)
			for i in range(nvar):
				start[i] = np.polyval(pol[:, i], prm[0])

		if loop_counter > 0:
			print('predictor : ', start, fit_type)
			for i in range(len(varia)):
				pyqcm.set_parameter(varia[i], start[i])
				print(' ---> ', varia[i], ' = ', start[i])
		try:
			pyqcm.banner('loop index = {:d}, {:s} = {:1.4f}'.format(loop_counter + 1, loop_param, prm[0]), '=', skip=1)
			func()
			P = pyqcm.parameters()
			for i in range(len(varia)):
				sol[0, i] = P[varia[i]]

		except pyqcm.OutOfBoundsError as E:
			print('variable ', E.variable + 1, ' is out of bounds: abs(', varia[E.variable], ') > ', max[E.variable])
			if loop_counter == 0 or not adjust:
				print('Out of bound on starting values, aborting')
				return
			else:
				retry = True
				continue
		except pyqcm.TooManyIterationsError as E:
			if loop_counter == 0 or not adjust:
				print('Cannot converge on starting values, aborting')
				return 
			else:
				retry = True
				continue

		if loop_counter > 2 and not retry and adjust:
			if(__predict_good(sol[0, :], start, 0.01) and step <  loop_range[2]):
				step *= 1.5
				print('readjusting step to ', step)
			elif(__predict_bad(sol[0, :], start, 0.2) and step > loop_range[2] / 10):
				step /= 1.5
				print('readjusting step to ', step)

		loop_counter += 1

		if control_func is not None:
			if not control_func():
				if loop_counter < 2 or not adjust:
					print('control failed on starting. aborting. Maybe try with different starting point')
					break
				retry = True
				continue
	if nretry > max_retry:
		pyqcm.banner('controlled loop ended on signal by control function', '%')
	else:
		pyqcm.banner('controlled loop ended normally', '%')

################################################################################
# performs a loop on a parameter while trying to keep a fixed density

def fixed_density_loop(
	mu,
	target_n,
	kappa=1.0,
	maxdmu=0.2,
	func=None,
	loop_param=None, 
	loop_values=None,
	var_param=None,
	dens_tol=0.001,
	dir='',
	measure=None,
	cluster_density=False
):
	"""Performs a loop while trying to keep a fixed density

	:param float mu: initial value of mu
	:param float target_n: desired value of the density
	:param kappa: initial guess of the compressibility
	:param maxdmu: maximum change in mu at each step (absolute value)
	:param func: function called at each step of the loop. No required argument. returns None.
	:param str loop_param: name of the parameter looped over
	:param float loop_values: an array of values of loop_param
	:param str var_param: array of variational parameters (names) to be predicted
	:param dens_tol: tolerance on the value of the density
	:param dir: directory of calling script
	:param measure: name of function to be called when the desired density is reached
	:param cluster_density: if True, uses the cluster density instead of the lattice density

	"""
	first_time = True
	print('fixed density loop over ', loop_param)
	if loop_param is None:
		raise pyqcm.MissingArgError('The name of the control parameter ("loop_param") is missing')
	if loop_values is None:
		raise pyqcm.MissingArgError('An array of values (loop_values) must be provided')

	mu_list = np.empty(len(loop_values))  # current solution + 3 previous solutions
	
	var = {}
	if var_param != None:
		for v in var_param:
			var[v] = np.empty(len(loop_values))

	i = 0
	restart = False
	if dir != '':
		restart = os.path.isfile(dir+'/density_loop.tsv')
	if restart:
		first_time = False
		dens_sol = np.genfromtxt(dir+'/density_loop.tsv', names=True)
		last_mu = dens_sol[loop_param][-1]
		i = 0
		for j, x in enumerate(dens_sol[loop_param]):
			if np.abs(x - last_mu) < 1e-8:
				i = j
		# i = loop_values.tolist().index(dens_sol[loop_param][-1])
		if i != len(dens_sol)-1:
			print('mismatch between the number of solutions in file ({:d}) and the index of loop parameter in loop_values ({:d})'.format(len(dens_sol),i+1))
		mu_list[0:i+1] = dens_sol['mu']

	while i < len(loop_values):
		pyqcm.set_parameter(loop_param, loop_values[i])
	
		# predicting the starting value from previous solutions (do nothing if loop_counter=0)
		fit_type = ''
		if i == 1:
			mu = mu_list[0]  # the starting point is the previous value in the loop
			pyqcm.set_parameter('mu', mu)
		if i == 2 :
			fit_type = ' (linear fit)'
			mu = poly_predictor(loop_values, mu_list, i, 1)
			if var_param != None:
				for v in var_param:
					pyqcm.set_parameter(v, poly_predictor(loop_values, var[v], i, 1))
		elif i > 2 :
			fit_type = ' (quadratic fit)'
			mu = poly_predictor(loop_values, mu_list, i, 2)
			if var_param != None:
				for v in var_param:
					pyqcm.set_parameter(v, poly_predictor(loop_values, var[v], i, 2))
			print('**** predictor: U =', loop_values[i-3:i], '  & mu = ', mu_list[i-3:i], ' --->  (U, mu) = (', loop_values[i], ',', mu, ')')

		print('*'*80)
		print(loop_param, ' = {:1.4f}'.format(loop_values[i]))

		# density search:
		max_trials = 12
		mu2 = np.zeros(max_trials)
		n = np.zeros(max_trials)
		var2 = {}
		if var_param != None:
			for v in var_param:
				var2[v] = np.empty(max_trials)

		mu_converged = False
		for j in range(max_trials+1):
			#+++ CALLING THE FUNCTION +++
			pyqcm.set_parameter('mu', mu)
			pyqcm.banner('mu = {:1.4f}'.format(mu), '+')
			func()
			if cluster_density:
				dens = pyqcm.cluster_averages()['mu'][0]
			else:
				dens = pyqcm.averages()['mu']
			n[j] = target_n-dens
			print('density = ', dens, '\t delta n = ', n[j])
			#++++++++++++++++++++++++++++
			if np.abs(n[j]) < dens_tol:
				mu_converged = True
				break
			mu2[j] = mu
			P = pyqcm.parameters()
			if var_param != None:
				for v in var_param:
					var2[v][j] = P[v]

			if j==0:
				dmu = n[j]/kappa
				if np.abs(dmu) > maxdmu:
					dmu *= maxdmu/np.abs(dmu)
				mu += dmu    
			elif j==1:
				mu_new = poly_predictor(n, mu2, j+1, 1)
				dmu = mu_new-mu
				if np.abs(dmu) > maxdmu:
					dmu *= maxdmu/np.abs(dmu)
				mu += dmu    
				if var_param != None:
					for v in var_param:
						pyqcm.set_parameter(v, poly_predictor(n, var2[v], j+1, 1))
			else:
				mu_new = poly_predictor(n, mu2, j+1, 2)
				dmu = mu_new-mu
				if np.abs(dmu) > maxdmu:
					dmu *= maxdmu/np.abs(dmu)
				mu += dmu    
				if var_param != None:
					for v in var_param:
						pyqcm.set_parameter(v, poly_predictor(n, var2[v], j+1, 2))

		if not mu_converged:
			print('failed to find the value of chemical potential. Aborting')
			exit(1)
		mu_list[i] = mu

		des, val = pyqcm.properties()
		f = open('density_loop.tsv', 'a')
		if first_time:
			f.write(des + '\n')
			first_time = False
		f.write(val + '\n')
		f.close()
		if measure != None:
			measure()

		P = pyqcm.parameters()
		if var_param != None:
			for v in var_param:
				var[v][i] = P[v]
		i += 1    



def poly_predictor(X,Y,i,n):
	x = X[i-n-1:i]  # the n+1 previous values of X
	y = Y[i-n-1:i]  # the n+1 previous values of Y
	pol = np.polyfit(x, y, n)
	return np.polyval(pol, X[i])


################################################################################
def fade(F, p1, p2, n):
	"""
	fades the model between two sets of parameters, in n steps

	:param F: task to perform wihtin the loop
	:param dict p1: first set of parameters
	:param dict p2: second set of parameters
	:param n: number of steps

	"""

	lambda_array = np.linspace(0.0, 1.0, n)
	for L in lambda_array:
		par = {}
		for p in p2:
			par[p] = L*p2[p]
			if p in p1:
				par[p] += (1-L)*p1[p]
		for p in p1:
			if p not in p2:
				par[p] = (1-L)*p1[p]
		
		for p in par:
			pyqcm.set_parameter(p, par[p])
		pyqcm.new_model_instance()
		F()


################################################################################
def Hartree(F, couplings, maxiter=10, eps_algo=0):
	"""Performs the Hartree approximation

	:param F: task to perform wihtin the loop
	:param [class hartree] couplings: list of Hartree couplings (or single coupling)
	:param int maxiter: maximum number of Hartree iterations
	:param int eps_algo: number of elements in the epsilon algorithm convergence accelerator = 2*eps_algo + 1 (0 = no acceleration)

	"""

	global first_time

	pyqcm.banner('Hartree procedure', c='*', skip=1)
	hartree_converged = False
	diff_tot = 1e6
	var_data = np.empty((len(couplings), maxiter+2))

	if eps_algo:
		for C in couplings:
			C.init_epsilon(maxiter, eps_algo)

	iter = 0
	while True:
		F()
		iter += 1
		pyqcm.banner('Hartree iteration {:d}'.format(iter), '-')
		diff_tot = 0
		hartree_converged = True
		for i,C in enumerate(couplings):
			C.update()
			C.print()
			diff_tot += np.abs(C.diff)
			hartree_converged = hartree_converged and C.converged()
			var_data[i, iter] = C.vm


		print('total difference = {:g}'.format(diff_tot))

		if hartree_converged:
			pyqcm.write_summary('hartree.tsv', first = first_time)
			first_time = False
			break

		if iter > maxiter :
			print('Maximum number of Hartree iterations exceeded! Aborting...')
			exit(1)

	pyqcm.banner('Hartree procedure has converged', c='*')
	


################################################################################
# performs a loop (a linear trajectory) between two subsets of parameters
def linear_loop(
	N,
	func=None,
	varia=None,
	params = None,
	predict=True
):
	"""Performs a controlled loop for VCA or CDMFT with a predictor

	:param N: number of intervals within the loop
	:param func: function called at each step of the loop
	:param [str] varia: names of the variational parameters
	:param {str:(float,float)} P: dict of parameters to vary with initial and final values 
	:param boolean predict: if True, uses a linear or quadratic predictor

	"""

	pyqcm.banner('linear loop over {:d} values'.format(N), '%', skip=1)
	if type(params) != dict:
		raise TypeError('the argument "P" of linear_loop must be a dict of 2-tuples')
	if varia is None:
		raise pyqcm.MissingArgError('variational parameters (varia) must be specified')
	if type(varia) != list:
		raise TypeError('the argument "varia" of controlled_loop must be a list')

	current_value = {}

	nvar = len(varia)  # number of variational parameters
	sol = np.empty((4, nvar))  # current solution + 3 previous solutions
	start = [0]*nvar
	for iter in range(N+1):
		# sets the control parameters of the loop
		for x in params:
			current_value[x] = params[x][0]*(N-iter)/N + params[x][1]*iter/N
			pyqcm.set_parameter(x, current_value[x])
		
		sol = np.roll(sol, 1, 0)  # cycles the solutions
		# predicting the starting value from previous solutions (do nothing if loop_counter=0)
		fit_type = ''
		if iter == 1 or not predict:
			start = sol[1, :]  # the starting point is the previous value in the loop
		elif iter > 1:
			if iter == 2:
				fit_type = ' (linear fit)'
				x = [iter-1, iter-2]  # the two previous values of the loop parameter
				y = sol[1:3, :]  # the two previous values of the loop parameter
				pol = np.polyfit(x, y, 1)
			else:
				fit_type = ' (quadratic fit)'
				x = [iter-1, iter-2, iter-3]  # the three previous values of the loop parameter
				y = sol[1:4, :]  # the three previous values of the loop parameter
				pol = np.polyfit(x, y, 2)
			for i in range(nvar):
				start[i] = np.polyval(pol[:, i], iter)

		pyqcm.banner('loop index = {:d}'.format(iter + 1), '=')
		for x in params:
			print(x, ' ===> {:1.4f}'.format(current_value[x]))
		if iter > 0:
			print('predictor: ', start, fit_type)
		try:
			func()
			P = pyqcm.parameters()
			for i in range(len(varia)):
				sol[0, i] = P[varia[i]]

		except pyqcm.OutOfBoundsError as E:
			print('variable ', E.variable + 1, ' is out of bounds: abs(', varia[E.variable], ') > ', max[E.variable])
			print('Out of bound on starting values, aborting')
			return
		except pyqcm.TooManyIterationsError as E:
			print('Cannot converge on starting values, aborting')
			return 

	pyqcm.banner('linear loop ended normally', '%')
