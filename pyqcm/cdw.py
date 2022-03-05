import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(precision=10, linewidth=512, suppress=True)

################################################################################
class cluster:
    def __init__(self, sites, e):
        """
        This class is used to define a supercluster and the associated superlattice
        :param [[int]] sites: list of sites of the cluster (integer components), or tuple of such collections if more than one cluster in a supercluster.
        :param [[int]] e: superlattice basis vectors (integer components, the number of vectors is the dimension). Must be right oriented.
        """
        e = np.array(e, dtype='int')
        self.e = e
        self.dim = self.e.shape[0]
        assert (self.dim > 0 and self.dim <= 3), 'the dimension must be 1, 2 or 3'
        self.E = np.eye(3, dtype='int')
        self.E[0:self.dim, :] = self.e
        S = self.E
        self.vol = int(np.rint(np.linalg.det(S)))  # fixed : go to nearest integer
        assert self.vol > 0, 'Please use a lattice with a positive orientation. Aborting'

        self.minors = np.empty((3,3), dtype=int)
        self.minors[0,0] =  S[1,1]*S[2,2] - S[2,1]*S[1,2]
        self.minors[1,0] = -S[1,0]*S[2,2] + S[2,0]*S[1,2]
        self.minors[2,0] =  S[1,0]*S[2,1] - S[1,1]*S[2,0]
        self.minors[0,1] = -S[0,1]*S[2,2] + S[2,1]*S[0,2]
        self.minors[1,1] =  S[0,0]*S[2,2] - S[2,0]*S[0,2]
        self.minors[2,1] = -S[0,0]*S[2,1] + S[0,1]*S[2,0]
        self.minors[0,2] =  S[0,1]*S[1,2] - S[1,1]*S[0,2]
        self.minors[1,2] = -S[0,0]*S[1,2] + S[1,0]*S[0,2]
        self.minors[2,2] =  S[0,0]*S[1,1] - S[0,1]*S[1,0]
        # print('minors:\n', self.minors)
        print('super unit cell volume = ', self.vol)

        # folding the sites
        if type(sites) is not tuple:
            sites = (sites,)
        sitesF = []
        sitesO = []
        clus=[]
        for i,c in enumerate(sites):
            csites = np.array(c, dtype='int')
            for x in csites:
                sitesF += [self.fold_shifted(x)[1]]
                sitesO += [x]
                clus += [i]
        self.sitesF = np.array(sitesF, dtype='int') # sites, folded into a conventional unit cell
        self.sites = np.array(sitesO, dtype='int') # sites, folded into a conventional unit cell
        self.clus = np.array(clus, dtype='int') # cluster associated with each site
        self.N = self.sites.shape[0]
        print('liste of sites and folded sites:')
        for i in range(self.N):
            print(i+1, '\t', self.sites[i], '\t', self.sitesF[i], '\tdiff = ', self.sitesF[i]-self.sites[i])
        # dictionary of sites (for getting the index from the position)
        self.siteF_index = {}
        for i,c in enumerate(sitesF):
            S = '[{:d},{:d},{:d}]'.format(c[0], c[1], c[2])
            self.siteF_index[S] = i




    def fold_shifted(self, r):
        """
        folds an integer vector back to the shifted unit cell of the superlattice
        This shifted unit cell is such that all of its vectors have positive or zero components in terms of superlattice vectors

        :param [int] r: integer component vector to be folded into the unit cell
        :return (int, [int], [int]): I, S, R where r = S + R and R is a superlattice vector, I is the index of the folded site

        """
        r = np.array(r, dtype=int)
        R = np.empty(3, dtype=int)
        S = np.empty(3, dtype=int)
        Q = np.empty(3, dtype=int)

        R = r@self.minors
        assert np.linalg.norm(R@self.E  - r*self.vol) < 1e-6, 'Folding error here : r = {}, R={}, check={}'.format(r,R, self.vol*r@np.linalg.inv(self.E))

        # now r = (R[0] e[0] + R[1] e[1] + R[2] e[2]) + (Q[0] e[0] + Q[1] e[1] + Q[2] e[2])/self.vol
        Q[0] = R[0] % self.vol
        R[0] //= self.vol
        if(Q[0] < 0):
            Q[0] += self.vol
            R[0] -= 1

        Q[1] = R[1] % self.vol
        R[1] //= self.vol
        if(Q[1] < 0):
            Q[1] += self.vol
            R[1] -= 1

        Q[2] = R[2] % self.vol
        R[2] //= self.vol
        if(Q[2] < 0):
            Q[2] += self.vol
            R[2] -= 1
	
        assert np.linalg.norm(R@self.E + Q@self.E//self.vol - r) < 1e-6, 'Folding error here : r = {}, R={}, Q={}'.format(r,R,Q)

        # S = self.E[0]*Q[0] + self.E[1]*Q[1] + self.E[2]*Q[2]
        S = np.dot(Q,self.E)
        S[0] //= self.vol
        S[1] //= self.vol
        S[2] //= self.vol
        Q = R
        # R = self.E[0]*Q[0] + self.E[1]*Q[1] + self.E[2]*Q[2]
        R = np.dot(Q,self.E)

  
        # then compensate between R and S to make the folding d-dimensional if d<3
        if(self.dim < 3):
            R -= self.E[2]*Q[2]
            S += self.E[2]*Q[2]
        if(self.dim  < 2):
            R -= self.E[1]*Q[1]
            S += self.E[1]*Q[1]

        try:
            Sstr = '[{:d},{:d},{:d}]'.format(S[0], S[1], S[2])
            I = self.siteF_index[Sstr]
        except:
            I = None

        assert np.linalg.norm(S+R-r)<1e-6, 'folding operation failed! r = {}, S={}, R={}'.format(r,S,R)
        return I, S, R

    def fold(self, r):
        """
        like 'fold_shifted', but this time folds back to the cluster itself, i.e., a different unit cell
        :param [int] r: integer component vector to be folded into the unit cell
        :return (int, [int], [int]): I, S, R where r = S + R and R is a superlattice vector, I is the index of the folded site
        """
        I, S, R = self.fold_shifted(r)
        if I != None:
            S += (self.sites[I] - self.sitesF[I])
            R += (self.sitesF[I] - self.sites[I])
        # print('fold : r = ', r, ' I = ', I, '  S = ', S, '  R = ', R)
        return I, S, R

    def draw(self, basis=None, plt_ax=None):
        """
        Plots the sites, the shifted sites and the superlattice vectors
        :param [[float]] basis: the real space, geometric basis

        """
        assert self.dim == 2, 'the draw() function only works in two dimensions'

        if plt_ax is None:
            ax = plt.gca()
        else:
            ax = plt_ax

        if basis is None:
            basis = np.eye(3)
        
        e = self.e
        e = e@basis
        ax.set_aspect(1)
        ax.plot([0, e[0,0]], [0, e[0,1]], 'r-', lw=0.5)
        ax.plot([0, e[1,0]], [0, e[1,1]], 'r-', lw=0.5)
        S = self.sites
        S = S@basis
        ax.plot(S[:,0], S[:,1], 'bo', ms=6)
        S = self.sitesF
        S = S@basis
        ax.plot(S[:,0], S[:,1], 'go', ms=3)
        if plt_ax is None:
            plt.show()

    def draw_cdw(self, cdw, basis=None, plt_ax=None):
        """
        Plots the amplitudes of a cdw
        :param [float] cdw: the cdw amplitudes (same order as the sites)
        :param [[float]] basis: the real space, geometric basis

        """
        assert self.dim == 2, 'the draw() function only works in two dimensions'

        if plt_ax is None:
            ax = plt.gca()
        else:
            ax = plt_ax

        if basis is None:
            basis = np.eye(3)
        
        e = self.e
        e = e@basis
        ax.set_aspect(1)
        ax.plot([0, e[0,0]], [0, e[0,1]], 'r-', lw=0.5)
        ax.plot([0, e[1,0]], [0, e[1,1]], 'r-', lw=0.5)
        S = self.sites
        S = S@basis
        fac = 12/np.max(np.abs(cdw))
        for j1 in range(-1,2,1):
            for j2 in range(-1,2,1):
                es = j1*e[0,:] + j2*e[1,:]
                for i in range(S.shape[0]):
                    ss = es + S[i,:]
                    if cdw[i] > 0: ax.plot(ss[0], ss[1], 'bo', ms=fac*cdw[i])
                    else: ax.plot(ss[0], ss[1], 'ro', ms=-fac*cdw[i])
        if plt_ax is None:
            plt.show()

################################################################################

def cdw_energy(C, U, _V, n, cluster=False, pr=False):
    """
    Computes the potential energy per site associated with a given CDW pattern

    :param (class cluster) C: periodic cluster
    :param U : on-site interaction
    :param [([int], float)] V: density-density interactions
    :param [float] n: density pattern
    :param boolean cluster: If True, limits the computation to the cluster (no inter-cluster interactions)
    :param boolean pr: If True, prints progress
    :return float: potential energy

    """
    
    V = []
    E = 0.0
    n = np.array(n)
    assert type(_V) == list
    for y in _V:
        assert type(y) == tuple

    for v in _V:
        V += [(np.array(v[0], dtype=int), v[1])]

    # contribution of the on-site energy
    for i in range(C.N):
        E += 0.25*U*n[i]*n[i]

    if pr:
        print('\nComputing the CDW energy with n = ', n)
    for v in V:
        if cluster:
            for i,x in enumerate(C.sites):
                R = x + v[0]
                I = None
                for j,r in enumerate(C.sites):
                    if np.linalg.norm(R-r) < 1e-6:
                        E += n[i]*n[j]*v[1]
                        break

        else:
            for i,x in enumerate(C.sitesF):
                j, S, R = C.fold(x + v[0])
                if j != None:
                    ener = n[i]*n[j]*v[1]
                    if np.abs(ener) > 1e-6 and pr:
                        print(v, '(', i+1, ',', j+1, ') --> ', ener)
                    E += ener
            
    return E/C.N


#----------------------------------------------------------------------
def cdw_eigenstates(C, _V, plt_ax=None, basis=np.eye(3)):
    """
    Computes the possible CDW states of the cluster, to be used with the Hartree approximation

    :param cluster C: periodic cluster
    :param [([int], float)] V: density-density interactions
    :return ([float], [float,float], [float,float]): the array of eigenvalues, the matrix of eigenvectors, and the inter-cluster interaction matrix

    """

    V = []
    E = 0.0
    Vic = np.zeros((C.N, C.N))
    Vc = np.zeros((C.N, C.N))
    S2 = C.sites@basis
    if plt_ax:
        plt_ax.set_aspect(1.0)
        plt_ax.scatter(S2[:,0], S2[:,1], [24]*S2.shape[0],color='gray')

    assert type(_V) == list
    for y in _V:
        assert type(y) == tuple

    for v in _V:
        V += [(np.array(v[0], dtype=int), v[1])]

    # print('\ninter-cluster V matrix contributions:\n')
    print('-'*80)
    for v in V:
        for i,x in enumerate(C.sites):  # loop over site 1
            j, S, R = C.fold(C.sites[i] + v[0])
            if j != None:
                if C.clus[i] == C.clus[j] and np.linalg.norm(R) < 1e-6:
                    Vc[i,j] += v[1]
                else:
                    Vic[i,j] += v[1]
                if plt_ax != None:
                    DX = v[0]@basis
                    plt_ax.plot([S2[i,0], S2[i,0]+DX[0]], [S2[i,1], S2[i,1]+DX[1]])
            j, S, R = C.fold(C.sites[i] - v[0])
            if j != None:
                if C.clus[i] == C.clus[j] and np.linalg.norm(R) < 1e-6:
                    Vc[i,j] += v[1]
                else:
                    Vic[i,j] += v[1]
                if plt_ax != None:
                    DX = v[0]@basis
                    plt_ax.plot([S2[i,0], S2[i,0]-DX[0]], [S2[i,1], S2[i,1]-DX[1]])

    print('intra-cluster V matrix:\n',Vc)
    print('inter-cluster V matrix:\n',Vic)

    w, v = np.linalg.eigh(Vic)
    w = np.round(w,10)
    for i in range(C.N):
        # y = v[:,i]
        # for j in range(C.N):
        #     if np.abs(y[j] > 1e-2):
        #         break
        # y /= y[j]
        # v[:,i] = y
        print('\neigenvalue ', w[i], ' :\n', v[:, i])
    return w, v, Vic, Vc



#----------------------------------------------------------------------
def draw_V(C, _V, plt_ax, basis=None):
    """
    draws the links associated with extended interactions

    :param (class cluster) C: periodic cluster
    :param [([int], float)] V: density-density interactions
    :param plt_ax: matplotlib axis set; triggers plotting
    :param [[float]] basis: the real space, geometric basis
    :return: None

    """
    
    V = []
    assert type(_V) == list
    for y in _V:
        assert type(y) == tuple

    for v in _V:
        V += [(np.array(v[0], dtype=int), v[1])]

    C.draw(basis, plt_ax)
    S2 = C.sites@basis
    color=['r', 'g', 'b', 'y', 'k']
    ic = 0
    for v in V:
        for i,x in enumerate(C.sitesF):
            j, S, R = C.fold_shifted(x + v[0])
            if j != None:
                plt_ax.plot([S2[i,0],S2[j,0]], [S2[i,1],S2[j,1]], '-', c = color[ic%5])
        ic += 1
    plt.show()

#----------------------------------------------------------------------
def draw_Vic(C, _V, plt_ax, basis=None):
    """
    draws the links associated with extended interactions

    :param (class cluster) C: periodic cluster
    :param [([int], float)] V: density-density interactions
    :param plt_ax: matplotlib axis set; triggers plotting
    :param [[float]] basis: the real space, geometric basis
    :return: None

    """
    
    V = []
    assert type(_V) == list
    for y in _V:
        assert type(y) == tuple

    for v in _V:
        V += [(np.array(v[0], dtype=int), v[1])]

    C.draw(basis, plt_ax)
    S2 = C.sites@basis
    color=['r', 'g', 'b', 'y', 'k']
    ic = 0
    for v in V:
        for i,x in enumerate(C.sitesF):
            j, S, R = C.fold_shifted(x + v[0])
            
            if j != None:
                plt_ax.plot([S2[i,0],S2[i,0]+DX[0]], [S2[i,1],S2[i,1]+DX[1]], '-', c = color[ic%5])
        ic += 1
    plt.show()

#----------------------------------------------------------------------
def draw_mode(C, X, plt_ax, basis=None):
    """
    draws the links associated with extended interactions

    :param (class cluster) C: periodic cluster
    :param [float] X: amplitude (one per site)
    :param plt_ax: matplotlib axis set; triggers plotting
    :param [[float]] basis: the real space, geometric basis
    :return: None

    """
    
    S2 = C.sites@basis
    plt_ax.scatter(S2[:,0], S2[:,1], [24]*S2.shape[0],color='gray')
    X = X/np.linalg.norm(X)
    X *= 48
    col = ['b']*len(X)
    for i in range(len(X)):
        if X[i] < 0 :
            col[i] = 'r'
    X = np.abs(X)
    plt_ax.scatter(S2[:,0], S2[:,1], X, c=col)

#----------------------------------------------------------------------
def sdw_energy(C, _J, sdw, cluster=False, pr=False):
    """
    Computes the energy per site associated with a given SDW pattern

    :param (class cluster) C: periodic cluster
    :param [([int], [[float]])] J: spin interactions
    :param [[float]] sdw: spin pattern (3 x N) : for each site, 3 angles are specifies : chi,theta,phi. The amplitude of the spin is sin(chi), its polar angle is theta and phi is the azimutal angle
    :param boolean cluster: If True, limits the computation to the cluster (no inter-cluster interactions)
    :param boolean pr: If True, prints progress
    :return float: energy

    """
    
    J = []
    E = 0.0
    assert type(_J) == list
    for y in _J:
        assert type(y) == tuple

    for v in _J:
        J += [(np.array(v[0], dtype=int), np.array(v[1]))]

    sdw = np.array(sdw)
    if pr:
        print('\nComputing the SDW energy with sdw = ', sdw)
    s = np.zeros((C.N,3))
    for i in range(C.N):
        s[i,0] = np.abs(np.sin(sdw[0,i]))*np.sin(sdw[1,i])*np.cos(sdw[2,i])
        s[i,1] = np.abs(np.sin(sdw[0,i]))*np.sin(sdw[1,i])*np.sin(sdw[2,i])
        s[i,2] = np.abs(np.sin(sdw[0,i]))*np.cos(sdw[1,i])
    for v in J:
        if cluster:
            for i,x in enumerate(C.sites):
                R = x + v[0]
                for j,r in enumerate(C.sites):
                    if np.linalg.norm(R-r) < 1e-6:
                        E += s[i]@v[1]@s[j]
                        break

        else:
            for i,x in enumerate(C.sitesF):
                j, S, R = C.fold(x + v[0])
                if j != None:
                    ener = s[i]@v[1]@s[j]
                    # if np.abs(ener) > 1e-6 and pr:
                    #     print(v, '(', i+1, ',', j+1, ') --> ', ener)
                    E += ener
            
    return E/C.N

def convert_to_spins(C, x, pr=False):
    s = np.zeros((C.N,3))
    sdw = np.reshape(x,(3,8))
    for i in range(C.N):
        s[i,0] = np.abs(np.sin(sdw[0,i]))*np.sin(sdw[1,i])*np.cos(sdw[2,i])
        s[i,1] = np.abs(np.sin(sdw[0,i]))*np.sin(sdw[1,i])*np.sin(sdw[2,i])
        s[i,2] = np.abs(np.sin(sdw[0,i]))*np.cos(sdw[1,i])

        if pr:
            print(C.sites[i], '\t', s[i,:])
    return s
