################################################################################
# This file contains functions implementing the Variational Cluster Approximation (VCA)
################################################################################

import numpy as np
import pyqcm
import time

#-------------------------------------------------------------------------------
# MPI matter
comm = None
root = False
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    if comm.Get_rank() == 0:
        root = True
        print('MPI available. Working on ', comm.Get_size(), ' concurrent processes')
except ImportError:
    root = True

first_time = True

################################################################################
# PRIVATE FUNCTIONS
################################################################################

################################################################################
# evaluation of many instances of a function in parallel

def __evalF(func, x):
    n = x.shape[0]
    y = np.empty(n)
    for i in range(n):
        y[i] = func(x[i,:])
    return y    

def __evalF_mpi(func, x):
    n = x.shape[0]
    nproc = comm.Get_size()
    rank = comm.Get_rank()
    y = np.zeros(n)
    y_tmp = np.zeros(n)
    k = 0  # index of points
    for i in range(n):
        if i % nproc == rank:
            y_tmp[i] = func(x[i,:])

    comm.Allreduce(y_tmp, y)
    return y


################################################################################
# quasi-Newton method function evaluation

def __QN_values(func, x, step):
    """Computes the function func of n variables at :math:`2n+1` points at x and at n pairs of points located at :math:`\pm` step from x in each direction. For use in the quasi-Newton method.

    :param func: function of n variables
    :param x: base point (n components NumPy array)
    :param float step: step to take in each direction in computing derivatives
    :returns: array of :math:`2n+1` values

    """
    n = len(x)
    N = 2 * n + 1
    xp = np.empty((N, n))
    for i in range(N):
        xp[i, :] = np.copy(x)
    k = 0
    for i in range(n):
        k += 1
        xp[k, i] += step[i]
        k += 1
        xp[k, i] -= step[i]

    if comm is None:
        return __evalF(func, xp)
    else:    
        return __evalF_mpi(func, xp)


################################################################################
# quasi-Newton method

def __quasi_newton(func=None, start=None, step=None, accur=None, max=10, gtol=1e-4, bfgs=False, max_iteration=30, max_iter_diff=None):
    """Performs the quasi_newton procedure
    
    :param func: a function of N variables
    :param [float] start: the starting values
    :param [float] step: the steps used to computed the numerical second derivatives
    :param [float] accur: the required accuracy for each variable
    :param [float] max: maximum absolute value of each parameter
    :param float gtol: the gradient tolerance (gradient must be smaller than gtol for convergence)
    :param boolean bfgs: True if the BFGS method is used, otherwise the symmetric rank-1 formula is used (default)
    :param int max_iterations: maximum number of iterations, beyond which an exception is raised
    :param float max_iter_diff: maximum step to make
    :return (float, [float], [[float]]): tuple of x (the solution), gradient (array, the value of the gradient), hessian (matrix, the Hessian matrix)

    """
    n = len(start)
    gradient = np.zeros(n)
    gradient0 = np.zeros(n)
    dx = np.zeros(n)
    y = np.zeros(n)
    ihessian = np.eye(n)  # inverse Hessian matrix
    iteration = 0
    x = start

    while iteration < max_iteration:
        x0 = x
        x, dx, gradient, ihessian = __quasi_newton_step(iteration, func, x0, step, gradient, dx, bfgs)
        iteration += 1


        if (np.linalg.norm(gradient) < gtol) :
            if root:
                print('convergence on gradient after ', iteration, ' iterations')
            break

        for i in range(n):
            if(np.abs(x[i]) > max[i]):
                raise pyqcm.OutOfBoundsError(variable=i, iteration=iteration)

        converged = True
        for i in range(n):
            if np.abs(dx[i]) > accur[i]:
                converged = False
                break
        if converged:
            if root:
                print('convergence on position after ', iteration, ' iterations')
            break

    if iteration == max_iteration:
        raise pyqcm.TooManyIterationsError(max_iteration)

    return x, gradient, ihessian

################################################################################
# quasi-Newton method (step)

def __quasi_newton_step(iteration = 0, func=None, x=None, step=None, gradient=None, dx=None, bfgs=False, max_diff=None):
    """Performs a step of the quasi_newton procedure
    
    :param int iter: iteration number
    :param func: a function of N variables
    :param [float] x: the starting values
    :param [float] step: the steps used to computed the numerical second derivatives
    :param [float] dx: the previous difference
    :param boolean bfgs: True if the BFGS method is used, otherwise the symmetric rank-1 formula is used (default)
    :param float max_diff: maximum step to make
    :return (float, float, [float], [[float]]): tuple of x (the new point), dx (the difference between the current and previous point), gradient (array, the value of the gradient), hessian (matrix, the Hessian matrix)

    """
    n = len(x)
    y = np.zeros(n)
    ihessian = np.eye(n)  # inverse Hessian matrix

    F = __QN_values(func, x, step)
    gradient0 = np.copy(gradient)
    for i in range(n):
        gradient[i] = (F[2 * i + 1] - F[2 * i + 2]) / (2 * step[i])

    if iteration == 0:
        for i in range(n):
            ihessian[i, i] = step[i] * step[i] / (F[2 * i + 1] + F[2 * i + 2] - 2 * F[0])

    else:
        y = gradient - gradient0
        if bfgs:  # BFGS method
            z1 = 1.0e-14 + np.dot(y, dx)
            z1 = 1.0 / z1
            X = np.eye(n)
            X -= z1 * np.kron(y, dx).reshape((n, n))
            ihessian = np.dot(X.T, np.dot(ihessian, X))
            ihessian += z1 * np.kron(dx, dx).reshape((n, n))
        else:  # symmetric rank 1
            u = np.dot(ihessian, y)  # $u = H y_k$
            dx -= u  # the $\Delta x_k - H_k y_k$ of Wikipedia
            z1 = 1.0e-14 + np.dot(y, dx)
            z1 = 1.0 / z1
            ihessian += z1 * np.kron(dx, dx).reshape((n, n))

    dx = -np.dot(ihessian, gradient)

    if max_diff is not None:
        norm_dx = np.linalg.norm(dx)
        if norm_dx > max_diff:
            ratio = max_diff / norm_dx
            dx *= ratio
            if root:
                s1 = 'Norm of QN iteration variation dx is out of bounds: {} > {}'.format(norm_dx, max_diff)
                s2 = 'Readjusting: dx = dx * {:.4}%'.format(ratio * 100)
                print(s1, '\n', s2)

    x += dx
    if root:
        print('QN iteration no ', iteration+1, '\t x = ', x)

    return x, dx, gradient, ihessian


################################################################################
# Newton-Raphson method :  computation of the Hessian

def __NR_Hessian(func, x, step):
    """Computes the function func of n variables at :math:`1 + 2n + n(n-1)/2` points in order to fit a quadratic form x and n pairs of points located at :math:`\pm` step from x in each direction
    
    :param func: the function
    :param [float] x: the base point 
    :param float step: step taken in each direction
    :return (float, [float], [[float]]): the value of the function, the gradient, and the Hessian

    """
    n = len(x)
    N = 1 + 2*n + (n*(n-1))//2
    y = np.zeros(N)
    hess = np.zeros((n,n))

    xp = np.zeros((N, n))
    for k in range(N):
        xp[k, :] = np.copy(x)
 
    k = 0
    # points located at +/- step
    for i in range(n):
        k += 1
        xp[k, i] += step[i]
        k += 1
        xp[k, i] -= step[i]

    # points located off diagonal
    for i in range(n):
        for j in range(i):
            k += 1
            xp[k, i] += 0.707*(2*(i%2)-1)*step[i]
            xp[k, j] += 0.707*(2*(j%2)-1)*step[j]

    if comm is None:
        y = __evalF(func, xp)
    else:    
        y = __evalF_mpi(func, xp)

    # setting fit coefficients
    A = np.empty((N,N))
    A[:, 1:n+1] = xp
    for k in range(N):
        A[k,0] = 1.0
        m = n+1
        for i in range(n):
            for j in range(i+1):
                A[k, m] = xp[k, i]*xp[k,j]
                if i == j :
                    A[k, m] = xp[k, i]*xp[k,j]*0.5
                else:
                    A[k, m] = xp[k, i]*xp[k,j]        
                m += 1

    # inverting
    A1 = np.linalg.inv(A)
    y = np.dot(A1, y)
    val = y[0]
    grad = y[1:n+1]
    val += np.sum(grad*x)
    m = n+1
    for i in range(n):
        for j in range(i+1):
            hess[i,j] = y[m]
            hess[j,i] = y[m]
            if i == j :
                val += x[i] * x[j] * y[m] * 0.5
            else:
                val += x[i] * x[j] * y[m]
            m += 1
    grad += np.dot(hess,x)
    return val, grad, hess


################################################################################
# Newton-Raphson method

def __newton_raphson(func=None, start=None, step=None, accur=None, max=10, gtol=1e-4, max_iteration=30, max_iter_diff=None):
    """Performs the Newton-Raphson procedure
    
    :param func: a function of N variables
    :param [float] start: the starting values
    :param [float] step: the steps used to computed the numerical second derivatives
    :param [float] accur: the required accuracy for each variable
    :param [float] max: maximum absolute value of each parameter
    :param float gtol: the gradient tolerance (gradient must be smaller than gtol for convergence)
    :param int max_iterations:  maximum number of iterations, beyond which an exception is raised
    :param float max_iter_diff: optional maximum value of the maximum step
    :returns (float, [float], [[float]]): the value of the function, the gradient, and the Hessian

    """
    n = len(start)
    gradient = np.zeros(n)
    dx = np.zeros(n)
    y = np.zeros(n)
    iteration = 0
    x = start
    step0 = step

    while iteration < max_iteration:
        iteration += 1

        gradient0 = np.copy(gradient)
        dx, F, gradient, hessian = __newton_raphson_step(func, x, step)


        if np.linalg.norm(gradient) < gtol:
            if root:
                print('convergence on gradient after ', iteration, ' iterations')
            break

        ihessian = np.linalg.inv(hessian)

        if max_iter_diff is not None:
            norm_dx = np.linalg.norm(dx)
            if norm_dx > max_iter_diff:
                ratio = max_iter_diff / norm_dx
                dx *= ratio
                s1 = 'Norm of NR iteration variation dx is out of bounds: {} > {}'.format(norm_dx, max_iter_diff)
                s2 = 'Readjusting: dx = dx * {:.4}%'.format(ratio * 100)
                if root:
                    print(s1, '\n', s2)
                    
        x += dx

        #redefining the steps
        step_multiplier = 2.0
        for i in range(n):
            step[i] = np.abs(dx[i])
            if step[i] > step0[i]:
                step[i] = 2.0*step0[i]
            if step[i] < step_multiplier*accur[i]:
                step[i] = step_multiplier*accur[i]

        if root:
            print('NR iteration no ', iteration, '\t x = ', x, '\t steps = ', step, '\t dx = ', dx)

        for i in range(n):
            if(np.abs(x[i]) > max[i]):
                raise pyqcm.OutOfBoundsError(variable=i, iteration=iteration)

        converged = True
        for i in range(n):
            if np.abs(dx[i]) > accur[i]:
                converged = False
                break
        if converged:
            if root:
                print('convergence on position after ', iteration, ' iterations')
            break

    if iteration == max_iteration:
        raise pyqcm.TooManyIterationsError(max_iteration)

    return x, gradient, ihessian

################################################################################
# Newton-Raphson step

def __newton_raphson_step(func=None, x=None, step=None):
    """Performs a Newton-Raphson step
    
    :param func: a function of N variables
    :param [float] x: the current value of the parameters (changed by the function)
    :param [float] step: the steps used to computed the numerical second derivatives
    :returns ([float], float, [float], [[float]]): the change in position, the value of the function, the gradient, and the Hessian

    """

    F, gradient, hessian = __NR_Hessian(func, x, step)
    ihessian = np.linalg.inv(hessian)
    dx = -np.dot(ihessian, gradient)

    return dx, F, gradient, ihessian


################################################################################
# PUBLIC FUNCTIONS
################################################################################
# performs the VCA

def vca(var2sef=None, names=None, start=None, steps=None, accur=None, max=None, accur_grad=1e-6, max_iter=30, max_iter_diff=None, NR=False, hartree=None):
    """Performs a VCA with the QN or NR method
    
    :param var2sef: function that converts variational parameters to model parameters
    :param [str] names: names of the variational parameters
    :param [float] start: starting values
    :param [float] steps: initial steps
    :param [float] accur: accuracy of parameters (also step for 2nd derivatives)
    :param [float] max: maximum values that are tolerated
    :param float accur_grad: max value of gradient for convergence
    :param int max_iter: maximum number of iterations in the procedure
    :param float max_iter_diff: optional maximum value of the maximum step in the quasi-Newton method
    :param boolean NR: True if the Newton-Raphson method is used, False if the quasi-Newton method is used
    :param (class hartree) hartree: Hartree approximation couplings (see pyqcm/hartree.py)
    :return: None
    
    """
    global first_time
    pyqcm.new_model_instance()
    L = pyqcm.model_size()[0]
    pyqcm.first_SEF = True

    if names is None:
        print('missing argument names : variational parameters must be specified')
        raise pyqcm.MissingArgError('names')

    nvar = len(names)  # number of variational parameters

    if max is None:
        max = 10*np.ones(nvar)
    elif len(max) != nvar:
        print('the argument "max" should have ', nvar, ' components')
        raise pyqcm.MissingArgError('max')

    if start is None:
        assert var2sef == None, 'the start argument is missing in vca(). Should be specified if var2sef is not None'
        start = [0.0]*nvar
        P = pyqcm.parameters()
        for i,v in enumerate(names):
            start[i] = P[v]

    elif len(start) != nvar:
        print('the argument "start" should have ', nvar, ' components')
        raise pyqcm.MissingArgError('start')

    if accur is None:
        accur = 1e-4*np.ones(nvar)

    if steps is None:
        steps = 10*np.array(accur)
    elif len(steps) != nvar:
        print('the argument "steps" should have ', nvar, ' components')
        raise pyqcm.MissingArgError('steps')

    def var2x(x):
        if var2sef is None:
            for i in range(len(names)): 
                pyqcm.set_parameter(names[i], x[i])
            print('x = ', x) # new
        else:
            var2sef(x)    
        pyqcm.new_model_instance()
        return pyqcm.Potthoff_functional(hartree)
        
    if hartree == None:
        pyqcm.banner('VCA procedure', '*')
    else:
        pyqcm.banner('VCA procedure (combined with Hartree procedure)', '*')
    var_val = pyqcm.__varia_table(names,start)
    print(var_val)


    try:
        if NR :
            sol, grad, iH = __newton_raphson(var2x, start, steps, accur, max, accur_grad, max_iteration=max_iter, max_iter_diff=max_iter_diff)  # Newton-Raphson process
        else:
            sol, grad, iH = __quasi_newton(var2x, start, steps, accur, max, accur_grad, False, max_iteration=max_iter, max_iter_diff=max_iter_diff)  # quasi-Newton process
    except pyqcm.OutOfBoundsError as E:
        print('variable ', E.variable + 1, ' is out of bounds: abs(', names[E.variable], ') > ', max[E.variable])
        raise pyqcm.OutOfBoundsError(variable=E.variable, iteration=E.iteration)

    except pyqcm.TooManyIterationsError as E:
        print('quasi-Newton method failed to converge after ', E.max_iteration, ' iterations')
        raise pyqcm.TooManyIterationsError(E.max_iteration)

    omega = var2x(sol)  # final, converged value
    if root:
        print('saddle point = ', sol)
        print('gradient = ', grad)
        print('second derivatives :', 1.0 / np.diag(iH))
        print('computing properties of converged solution...')
        print('omega = ', omega)
    ave = pyqcm.averages()


    # writes the solution in the standard file
    if root:
        H = np.linalg.inv(iH)  # Hessian at the solution (inverse of iH)
        val = ''
        for i in range(nvar):
            val += str(H[i, i]) + '\t'
        des = ''
        for i in range(nvar):
            des += '2der_' + names[i] + '\t'
        pyqcm.write_summary('vca.tsv', first = first_time, suppl_descr = des, suppl_values = val)
        first_time = False

    if root:
        pyqcm.banner('VCA ended normally', '*')

    return sol, 1.0/np.diag(iH)

################################################################################
def plot_sef(param, prm, accur_SEF=1e-4, hartree=None, show=True):
    """Draws a plot of the Potthoff functional as a function of a parameter param taken from the list prm. The results are going to be appended to 'sef.tsv'
    
    :param str param: name of the parameter (independent variable)
    :param [float] prm: list of values of the parameter
    :param float accur_SEF: precision of the computation of the self-energy functional
    :param (class hartree) hartree: Hartree approximation couplings (see pyqcm/hartree.py)
    :param boolean show: if True, the plot is shown on the screen.
    :returns: None

    """
    L = pyqcm.model_size()[0]
    pyqcm.first_SEF = True

    pyqcm.set_global_parameter('accur_SEF', accur_SEF)
    omega = np.empty(len(prm))
    for i in range(len(prm)):
        pyqcm.set_parameter(param, prm[i])
        pyqcm.new_model_instance()
        # print(pyqcm.parameter_set(opt='report'))
        # print('.'*80, '\n', pyqcm.cluster_parameters())
        omega[i] = pyqcm.Potthoff_functional(hartree)

        print("omega(", prm[i], ") = ", omega[i])
    
    if show:
        import matplotlib.pyplot as plt
        plt.xlim(prm[0], prm[-1])
        plt.plot(prm, omega, 'b-')
        plt.xlabel(param)
        plt.ylabel(r'$\omega$')
        plt.axhline(omega[0], c='r', ls='solid', lw=0.5)
        plt.title(pyqcm.parameter_string())
        plt.show()




################################################################################
def plot_GS_energy(param, prm, clus=0, file=None, plt_ax=None, **kwargs):
    """Draws a plot of the ground state energy as a function of a parameter param taken from the list `prm`. The results are going to be appended to 'GS.tsv'
    
    :param str param: name of the parameter (independent variable)
    :param [float] prm: list of values of the parameter
    :param int clus: label of the cluster (starts at 0)
    :param str file: if not None, saves the plot in a file with that name
    :param plt_ax: optional matplotlib axis set, to be passed when one wants to collect a subplot of a larger set
    :param kwargs: keyword arguments passed to the matplotlib 'plot' function
    :returns: None

    """
    import matplotlib.pyplot as plt
    if plt_ax == None:
        plt.figure()
        plt.gcf().set_size_inches(13.5/2.54, 9/2.54)
        ax = plt.gca()
    else:
        ax = plt_ax

    omega = np.empty(len(prm))
    for i in range(len(prm)):
        pyqcm.set_parameter(param, prm[i])
        pyqcm.new_model_instance()
        omega[i] = pyqcm.ground_state()[clus][0]
        print("omega(", prm[i], ") = ", omega[i])

        # writing the parameters in a progress file
        f = open('GS.tsv', 'a')
        des, val = pyqcm.properties()
        if i == 0:
            f.write('\n\n')
            f.write(des + '\n')
        f.write(val + '\n')
        f.close()

    
    ax.axhline(omega[0], c='r', ls='solid', lw=0.5)
    ax.plot(prm, omega, 'b-')
    ax.set_xlim(prm[0], prm[-1])
    if plt_ax == None:
        ax.set_xlabel(param)
        ax.set_ylabel('GS energy')
        ax.set_title(pyqcm.parameter_string())

    if file is not None:
        plt.savefig(file)
        plt.close()
    elif plt_ax == None:
        plt.show()


################################################################################
# performs the VCA

def vca_min(names=None, start=None, steps=None, accur=1e-4, ftol=1e-8, method='Nelder-Mead', hartree=None):
    """Performs the VCA assuming that the solution is a minimum of the Potthoff functional
    Uses minimization routines from scipy.optimize.
    
    :param [str] names: names of the variational parameters
    :param [float] start: starting values 
    :param [float] steps: initial steps (relevant to some minimization methods)
    :param float accur: accuracy of parameters
    :param float ftol: convergence criterion for the value of the SEF
    :param str method: minimization method used in scipy.optimize.minimize()
    :param (class hartree) hartree: Hartree approximation couplings (see pyqcm/hartree.py)
    :return: None

    """


    if names is None:
        print('missing argument names : variational parameters must be specified')
        raise pyqcm.MissingArgError('names')

    nvar = len(names)  # number of variational parameters
    L = pyqcm.model_size()[0]

    if start is None:
        print('missing argument start')
        raise pyqcm.MissingArgError('start')
    elif len(start) != nvar:
        print('the argument "start" should have ', nvar, ' components')
        raise pyqcm.MissingArgError('start')

    if accur is None:
        accur = 1e-4*np.ones(nvar)

    if steps is None:
        steps = 10*np.array(accur)
    elif len(steps) != nvar:
        print('the argument "steps" should have ', nvar, ' components')
        raise pyqcm.MissingArgError('steps')

    def F(x):
        for i in range(len(names)): 
            pyqcm.set_parameter(names[i], x[i])
            print('x = ', x) # new
        pyqcm.new_model_instance()
        return pyqcm.Potthoff_functional(hartree)

    from scipy.optimize import minimize
    ftol = 1e-4
    X = None

    if method == 'Nelder-Mead':
        initial_simplex = np.zeros((nvar+1,nvar))
        for i in range(nvar+1):
            initial_simplex[i, :] = start
        for i in range(nvar):
            initial_simplex[i+1, i] += steps[i]
        sol = minimize(F, start, method='Nelder-Mead', options={'maxfev':200, 'xatol': accur, 'fatol':ftol, 'initial_simplex': initial_simplex, 'adaptive': True, 'disp':True})
        iter_done = sol.nit
        X = sol.x

    elif method == 'Powell':
        sol = minimize(F, start, method='Powell', tol = ftol, options={'xtol': accur, 'ftol': ftol, 'disp':True})
        iter_done = sol.nit
        X = sol.x
        if type(X) != list:
            X = [X]

    elif method == 'CG':
        sol = minimize(F, start, method='CG', jac=False, tol = ftol, options={'eps':accur, 'maxiter':200, 'disp':True})
        iter_done = sol.nit
        X = sol.x

    elif method == 'BFGS':
        sol = minimize(F, start, method='BFGS', jac=False, tol = ftol, options={'eps':accur, 'disp':True})
        iter_done = sol.nit
        X = sol.x

    elif method == 'COBYLA':
        sol = minimize(F, start, method='COBYLA', options={'rhobeg':steps[0], 'maxiter':200, 'tol': ftol, 'disp':True})
        iter_done = sol.nfev
        X = sol.x
        
    else:
        print('unknown method specified for minimization: ', method)
        exit()   

    if sol.success == False:
        print(sol.message)
        raise pyqcm.MinimizationError()

    omega = F(X)  # final, converged value
    ave = pyqcm.averages()
    des, val = pyqcm.properties()

    val += time.strftime("%Y-%m-%d@%H:%M", time.localtime())
    # writes the solution in the standard file
    if root:
        fout = open('vca.tsv', 'a')
        des += 'time'    
        fout.write(des + '\n')
        fout.write(val + '\n')
        fout.close()

    if root:
        print('*************** VCA ended normally ***************')


################################################################################
# detects a continuous phase transition

def __transition(varia, P, bracket, step=0.001, verb=False):
    """Detects a transition as a function of external parameter param by looking at the equality between 
    :math:`\Omega(h=s)` and :math:`\Omega(h=0)` where *h* is a single variational parameter (Weiss field)
    and *s* is a step. 
    
    :param str varia: name of the variational parameter
    :param str P: name of the control parameter (the parameter that controls the transition)
    :param (float,float) bracket: bracketing values for parameter P that enclose the transition
    :param float step: small, but finite value *s* of the Weiss field 
    :param boolean verb: If True, prints progress
    :returns: the value of P at the transition

    """

    from scipy.optimize import brentq

    def F(x):
        pyqcm.set_parameter(P, x)
        pyqcm.set_parameter(varia, step)
        pyqcm.new_model_instance()
        Om1 = pyqcm.Potthoff_functional()
        pyqcm.set_parameter(varia, 1e-8)
        pyqcm.new_model_instance()
        Om0 = pyqcm.Potthoff_functional()
        if verb:
            print(P, '= ', x, '\tdelta Omega = ', Om1-Om0)
        return Om1-Om0


    x0, r = brentq(F, bracket[0], bracket[1], maxiter=100, full_output=True, disp=True)
    if r.converged == False:
        print('The root finding routine could not find a solution!')
        exit(1)
        
    return x0

################################################################################
# detects a continuous phase transition (meta function with loop over control parameter)

def transition_line(varia, P1, P1_range, P2, P2_range, delta, verb=False):
    """Builds the second-order transition line as a function of a control parameter P1. The results are written in the file `transition.tsv`

    :param str varia: variational parameter
    :param str P1: control parameter
    :param [float] P1_range: an array of values of P1
    :param str P2: dependent parameter
    :param (float) P2_range: 2-uple of values of P2 that bracket the transition
    :param boolean verb: If True, prints progress
    :param float delta: at each step, the new bracket for P2 will be P2c +/- delta, P2c being the previous critical value 
    :return: None

    """

    assert type(varia) == str, 'the first argument of transition_line must be the name of a single variational parameter'
    assert type(P1) == str, 'the second argument of transition_line must be the name of a single control parameter'
    assert type(P2_range) == list, 'the argument P2_range of transition_line must be a list of 2 floats'

    trans = np.empty((len(P1_range), 2))
    for i, p1 in enumerate(P1_range):
        pyqcm.set_parameter(P1, p1)
        p2c = __transition(varia, P2, P2_range, step=0.001, verb=verb)
        trans[i,0] = p1
        trans[i,1] = p2c
        P2_range[0] = p2c - delta
        P2_range[1] = p2c + delta
        print('---> transition at ', P1, ' = ', p1, '\t', P2, ' = ', p2c)
        np.savetxt('transition.tsv', trans[0:i+1, :], delimiter='\t', header=P1+'\t'+P2)
