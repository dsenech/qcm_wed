import numpy as np
import pyqcm

#-------------------------------------------------------------------------------
# global variables

var = []
w = []
wr = []
weight = []
first_time = True
first_time2 = True
Hyb = []
Hyb_down = []
hybrid_to_param = None

Gdim = 0
nmixed = 1
mixing = 0
clusters = None
maxfev = 500000

######################################################################
class DCA_cluster:
    #----------------------------------------------------------------------
    def __init__(self, name, sites, superlattice, lattice=None, basis=None):
        """
        :param str name: name of the cluster
        :param [[int]] sites: list of sites of the cluster, integer components
        :param [[int]] superlattice: superlattice basis vectors (integer components, the number of vectors is the dimension)
        :param [[float]]lattice: lattice basis vectors (integer components, optional). 
        :param [[float]] basis: basis vectors (real components). Optional. 

        """
        self.name = name
        self.sites = np.array(sites, dtype='int')
        self.n = self.sites.shape[0]  # number of sites
        self.superlattice = np.array(superlattice, dtype='int')
        self.dim = self.superlattice.shape[0]
        assert (self.dim > 0 and self.dim <= 3), 'DCA : the dimension must be 1, 2 or 3'

        if lattice == None:  # default value for lattice
            self.lattice = np.eye(3, dtype='int')[0:self.dim,:]
        else:
            self.lattice = np.array(lattice)
            assert(self.lattice.shape[1] == 3 and self.lattice.shape[0] == self.dim), 'The lattice vectors have the wrong dimension or number'

        if basis == None:   # default value for basis
            self.basis = np.eye(3)
        else:
            self.basis = np.array(basis)
            assert(self.basis.shape[1] == 3 and self.basis.shape[0] == 3), 'The physical basis or its vectors have the wrong dimension'
            assert np.dot(self.basis[0],self.basis[1]) < -1e-9, 'Please choose a basis with an obtuse angle'

        self.ibasis = np.linalg.inv(self.basis)
        self.metric = self.ibasis.T@self.ibasis
        print('metric in k-space:\n', self.metric)
        # S : superlattice vectors in terms of the lattice vectors
        S = np.eye(3, dtype='int')
        S[0:self.dim, :] = self.superlattice
        latt = np.eye(3)
        latt[0:self.dim, :] = self.lattice
        ilatt = np.linalg.inv(latt)
        print('lattice :\n', latt, '\n', ilatt)
        S = np.dot(S,ilatt.T)

        # Sd : basis vectors of the reciprocal superlattice
        self.vol = int(np.rint(np.linalg.det(S)))  # number of q-vectors
        print('volume of the supercell (in terms of the unit cell) = number of q-vectors :', self.vol)
        self.dual = np.eye(3, dtype='int')
        self.dual[0,:] = np.cross(S[1,:], S[2,:])
        self.dual[1,:] = np.cross(S[2,:], S[0,:])
        self.dual[2,:] = np.cross(S[0,:], S[1,:])
        print('basis vector in k-space:\n', self.dual[0:self.dim]/self.vol)

        # Q = wavevectors of the reciprocal superlattice within the original Brillouin zone (2D for the moment)
        Q = []
        for i in range(-self.vol, self.vol):
            for j in range(-self.vol, self.vol):
                for k in range(-self.vol, self.vol):
                    q = i*self.dual[0] + j*self.dual[1] + k*self.dual[2]
                    if (q[0] >= self.vol) or (q[1] >= self.vol) or (q[2] >= self.vol) or (q[0] < 0)  or (q[1] < 0)  or (q[2] < 0):
                        continue
                    Q.append(q)
        print('construction of ', len(Q), " wavevectors:")
        self.Q = np.array(Q, dtype='float')
        self.Q /= self.vol
        self.nQ = len(self.Q) # number of patches
        print("wavevectors of the reciprocal superlattice within the original Brillouin zone (x 2pi):\n", self.Q)

        q1 = np.array([1,0,0])
        q2 = np.array([0,1,0])
        q3 = np.array([0,0,1])

        if self.dim == 1:
            self.extension = [np.zeros(3), q1, -q1]
        elif self.dim == 2:
            self.extension = [np.zeros(3), q1, -q1, q2, -q2, q1+q2, -q1-q2, q1-q2, -q1+q2]
        elif self.dim == 3:
            self.extension = [np.zeros(3), q1, -q1, q2, -q2, q3, -q3, q1+q2, -q1-q2, q1-q2, -q1+q2,
            q1+q3, -q1-q3, q1-q3, -q1+q3, q3+q2, -q3-q2, q3-q2, -q3+q2, q1+q2+q3, -q1+q2+q3, q1-q2+q3, -q1-q2+q3, q1+q2-q3, -q1+q2-q3, q1-q2-q3, -q1-q2-q3]

        # creating the cluster models

        self.__translations()

    #----------------------------------------------------------------------
    def __fold(self, x):
        """
        Folds a position x (integer components) into the super-unit-cell. Returns the index of the site.
        """
        for i in range(self.n):
            diff = x - self.sites[i]
            xc = self.dual@diff
            xc, f = np.divmod(xc, self.vol)
            if np.linalg.norm(f) == 0:
                return i
        return self.vol

    #----------------------------------------------------------------------
    def __translations(self):
        """
        Analyzes the translations and builds the correspondence between wavevectors and symmetry representations
        to be built:
            self.irrep : the correspondence between irreducible representation and q-vector: self.irrep[i] is the 
            irrep number associated with q-vector no i
            self.gen : the symmetry generators, to be passed to pyqcm.new_cluster_model(...)
        """

        # constructing the generators

        self.irrep = np.array(range(self.vol))
        T = self.lattice[0]
        g = np.zeros(self.n, dtype='int')
        for i in range(self.n):
            g[i] = self.__fold(self.sites[i]+T)
        self.gen = [g+1]
        
        if self.__permutation_length(g) != self.vol:
            # now assume the second translation vector to be [0,1,0]
            T = self.lattice[1]
            g = np.zeros(self.n, dtype='int')
            for i in range(self.n):
                g[i] = self.__fold(self.sites[i]+T)
            self.gen += [g+1]

            if self.__permutation_length(self.gen[0]-1) + self.__permutation_length(g) != self.vol:
                T = self.lattice[2]
                g = np.zeros(self.n, dtype='int')
                for i in range(self.n):
                    g[i] = self.__fold(self.sites[i]+T)
                self.gen += [g+1]

        print('generators :', self.gen)

    #----------------------------------------------------------------------
    def __permutation_length(self, p):
        
        j = 0
        m = 1
        while p[j] != 0:
            j = p[j]
            m += 1
        return m

    #----------------------------------------------------------------------
    def draw_patches(self, center=True, lim=1.05):
        """
        :param boolean center: If True, centers the plot at (0,0)
        :param float lim: the plot goes from -lim to lim in each direction (lim = 1 is pi)
        """
        import matplotlib.pyplot as plt
        from scipy.spatial import Voronoi, voronoi_plot_2d
        # zone boundary
        if self.dim != 2:
            print('draw_patches works only in dimension 2')
            return
        ZB = np.array([[0,0],[2,0],[2,2],[0,2],[0,0]])

        ibasis = self.ibasis[0:2, 0:2]
        sites = np.dot(self.sites[:, 0:2],self.basis[0:2, 0:2])
        sup = np.dot(self.superlattice[0:2, 0:2],self.basis[0:2, 0:2])
        Qr = 2*np.dot(self.Q[:, 0:2], ibasis.T)
        ZB = np.dot(ZB,ibasis.T)

        # basis vectors of the reciprocal lattice
        latt = np.eye(3)
        latt[0:self.dim, :] = self.lattice
        evol = int(np.rint(np.linalg.det(latt)))
        print('volume of the unit cell:', evol)
        K = np.eye(3) # basis vectors of the reciprocal lattice
        K[0,:] = np.cross(latt[1,:], latt[2,:])
        K[1,:] = np.cross(latt[2,:], latt[0,:])
        K[2,:] = np.cross(latt[0,:], latt[1,:])
        K *= 2/evol
        K = K[:self.dim,:self.dim]
        print('basis of reciprocal space:')
        for i in range(self.dim):
            print(K[i])
        K = np.dot(K,ibasis.T)

        if center:
            a = 2
            F = np.array([[a,0],[0,a],[-a,a],[-a,0],[0,-a],[a,-a],[a,0]])
            F = np.dot(F,ibasis.T)
            Fc = np.zeros(F.shape)
            norms = np.array([0.25*np.dot(F[i],F[i]) for i in range(F.shape[0])])
            for i in range(F.shape[0]-1):
                Fc[i] = np.linalg.solve(0.5*F[[i,i+1],:], norms[i:i+2])
            Fc[F.shape[0]-1] = Fc[0]
            ZB = Fc
            print(ZB)
            for i in range(self.vol):
                q = Qr[i]
                for p in F:
                    if np.linalg.norm(q-p) < np.linalg.norm(q):
                        Qr[i,:] = q-p

        # plot in direct space
        ms = 8
        plt.subplot(1,2,1)
        plt.gca().set_aspect(1)
        plt.grid()
        plt.xlabel('$x$')
        plt.ylabel('$y$')
        plt.plot(sites[:,0], sites[:,1], 'ko', ms=ms)
        plt.arrow(0, 0, sup[0,0], sup[0,1], width=0.01, head_width=0.1, length_includes_head=True, color='r')
        plt.arrow(0, 0, sup[1,0], sup[1,1], width=0.01, head_width=0.1, length_includes_head=True, color='r')
        plt.plot([sup[0,0]], [sup[0,1]], color='w')
        plt.plot([sup[1,0]], [sup[1,1]], color='w')


        # graphique dans l'espace réciproque
        plt.subplot(1,2,2)
        plt.gca().set_aspect(1)

        plt.plot(Qr[:,0], Qr[:,1], 'ko', ms=ms)
        plt.grid()
        plt.xlabel('$k_x/\pi$')
        plt.ylabel('$k_y/\pi$')

        # graphique de Voronoi
        # extension de la liste des points sur les zones voisines
        QrV = []
        for i,q in enumerate(Qr):
            QrV.append(q)
            QrV.append(q + K[0])
            QrV.append(q + K[1])
            QrV.append(q - K[0])
            QrV.append(q - K[1])
            QrV.append(q - K[0] + K[1])
            QrV.append(q + K[0] - K[1])
            QrV.append(q - K[0] - K[1])
            QrV.append(q + K[0] + K[1])
        QrV = np.array(QrV)

        vor = Voronoi(QrV)
        voronoi_plot_2d(vor, plt.gca())

        if center:
            plt.xlim(-lim,lim)
            plt.ylim(-lim,lim)
        plt.plot(ZB[:,0],ZB[:,1], 'b-')
        plt.tight_layout()
        plt.show()

    #----------------------------------------------------------------------
    def distance(self, k, q):
        """ 
        computes the periodic distance betweem two wavevectors k and q, each normalized
        to be contained in the unit cube
        :param [float] k: first wavevector
        :param [float] q: second wavevector
        :returns float:

        """
        d = np.array([(k-q+qq)@self.metric@(k-q+qq) for qq in self.extension])
        return np.min(d)

        # S = np.sin((k-q)*np.pi)
        # # S = 1-np.cos((k-q)*2*np.pi)
        # return np.sqrt(S@self.metric@S)

    #----------------------------------------------------------------------
    def grid(self, N=16, grid_type='nearest', xi=0.05):
        """ 
        constructs the DCA momentum grid
        :param int N: number of wavevectors in each direction of the grid
        :param str grid_type: type of grid; one of 'nearest', 'smooth'
        :returns None:

        """

        # grid wavevectors are represented by a 1D array of real-valued vectors in the unit cube
        # in the list, the x component varies fastest, followed by the y-component, etc.
        self.Nk1D = N
        self.Nk = N
        if self.dim==2:
            self.Nk *= N 
        elif self.dim==3:
            self.Nk *= N*N
        self.k = np.zeros((self.Nk, 3))

        # I = x + N*(y + N*z), therefore  x = I%N, y = (I//N)%N,  z = I//(N*N)
        I = np.array(range(self.Nk))
        y, x = np.divmod(I, N)
        z, y = np.divmod(y, N)
        self.k[:,0] = x/N
        self.k[:,1] = y/N
        self.k[:,2] = z/N


        if self.dim == 1:
            self.k += np.array([1,0,0])/(2*N) # shift
        elif self.dim == 2:
            self.k += np.array([1,1,0])/(2*N) # shift
        else:
            self.k += np.array([1,1,1])/(2*N) # shift

        Q = self.Q

        self.weight = np.zeros((self.Nk,self.vol)) # weights of each grid point wrt to each cluster wavevector

        for i in range(len(self.k)):
            for j in range(self.nQ):
                self.weight[i,j] = self.distance(self.k[i,:], Q[j,:])
        self.weight = np.round(self.weight,6)
        # computing the periodic distance between each grid point and the cluster wavevectors
        if grid_type=='nearest':
            for i in range(self.Nk):
                j = np.argmin(self.weight[i,:])
                self.weight[i,:] = np.zeros(self.vol)
                self.weight[i,j] = 1.0
        elif grid_type=='smooth':
            dist = self.vol**(-2/self.dim)
            self.weight = np.exp(-self.weight * self.weight/(xi*dist))
        else:
            print('the grid type "'+grid_type+'" does not exist!')
            exit(1)
        # normalization
        for i in range(self.vol):
            self.weight[:,i] /= np.sum(self.weight[:,i])
        # print('weights of the grid points:\n', self.weight)

    #----------------------------------------------------------------------
    def plot_grid(self):
        import matplotlib.pyplot as plt
        if self.dim != 2:
            return
        ibasis = self.ibasis[0:2,0:2]
        Qr = 2*np.dot(self.Q[:, 0:2], ibasis.T)
        kr = 2*np.dot(self.k[:, 0:2], ibasis.T)
        ZB = 2*np.array([[0,0],[1,0],[1,1],[0,1],[0,0]])
        ZB = np.dot(ZB,ibasis.T)
        for i,q in enumerate(Qr):
            m = np.max(self.weight[:,i])
            rgb = np.zeros((self.Nk, 3))
            rgb[:,2] = self.weight[:,i]/m
            plt.gca().set_aspect(1)
            plt.plot(ZB[:,0],ZB[:,1], 'b-')
            plt.plot([q[0]], [q[1]], 'ro', ms=6)
            # plt.scatter(kr[:,0], kr[:,1], c=rgb)
            plt.scatter(kr[:,0], kr[:,1], s=16*(self.weight[:,i]/m))
            plt.title('q-vector :'+str(q))
            plt.show()
        plt.gca().set_aspect(1)
        W = np.sum(self.weight, axis=1)
        m = np.max(W)
        plt.plot(Qr[:,0], Qr[:,1], 'ro', ms=6)
        plt.scatter(kr[:,0], kr[:,1], s=16*W/m)
        plt.title('sum')
        plt.show()


    #----------------------------------------------------------------------
    def set_up_models(self, nb):

        pyqcm.set_global_parameter('periodic')
        no = self.n + nb
        self.varia = []
        for ik in range(self.vol):
            gen = np.zeros((len(self.gen), self.n + nb), dtype='int')
            for i in range(gen.shape[0]):
                gen[i, 0:self.n] = self.gen[i]
                gen[i, self.n:self.n+nb] = ik
            clus_name = 'clus{:d}'.format(ik)
            pyqcm.new_cluster_model(clus_name, self.n, nb, generators=gen, bath_irrep=True)
            for i in range(nb):
                pyqcm.new_cluster_operator(clus_name, 'eb{:d}'.format(i+1), 'one-body', [(self.n+i+1, self.n+i+1, 1), (self.n+i+1+no, self.n+i+1+no, 1)])
                self.varia += ['eb{:d}_1'.format(i+1)]
                e = []
                for j in range(self.n):
                    phase = self.Q[ik]@self.sites[j]*2*np.pi
                    z = np.exp(phase*1j)
                    e += [(j+1, self.n+i+1, z), (j+1+no, self.n+i+1+no, z)]
                pyqcm.new_cluster_operator_complex(clus_name, 'tb{:d}'.format(i+1), 'one-body', e)
                self.varia += ['tb{:d}_1'.format(i+1)]

    #----------------------------------------------------------------------
    def cluster_Green_function(self, z):
        for q in range(self.nQ):
            pass



###########################################################################
def __frequency_grid(type='sharp', beta='50', wc=2):
    """
    constructs a grid of frequencies along the imaginary axis for the distance function

    :param str type: type of grid ('sharp', 'ifreq', 'self')
    :param float beta : inverse temperature
    :param float wc : cutoff
    :returns str dist_function

    """
    global w, wr, weight

    wr = np.arange((np.pi / beta), wc + 1e-6, 2 * np.pi / beta)
    w = np.ones(len(wr), dtype=np.complex128)
    w = w * 1j
    w *= wr
    nw = len(w)
    if type == 'sharp':
        weight = np.ones(nw)
        weight *= 1.0 / nw
        dist_function = 'sharp_wc_{0:.1f}_b_{1:d}'.format(wc, int(beta))
    elif type == 'ifreq':
        weight = 1.0/wr
        weight *= 1.0 / weight.sum()
        dist_function = 'ifreq_wc_{0:.1f}_b_{1:d}'.format(wc, int(beta))
    elif type == 'self':
        weight = np.zeros(nw)
        Sig_inf = pyqcm.cluster_self_energy(0, 1.0e6j)
        for i, x in enumerate(w):
            Sig = pyqcm.cluster_self_energy(0, x) - Sig_inf
            weight[i] = np.linalg.norm(Sig)
        weight *= 1.0 / weight.sum()
        dist_function = 'self_wc_{0:.1f}_b_{1:d}'.format(wc, int(beta))
    else:
        raise pyqcm.WrongArgumentError('frequency_grid', type)
    return dist_function

######################################################################
def __set_Hyb(spin_down=False):
    """Computes the hybridization function

    :param boolean spin_down: if True, considers the spin down sector (when mixing=4)
    :returns: an array of arrays of matrices. Hyb[i], for cluster #i, is a (nw,d,d) Numpy array. with nw frequencies, and d sites.

    """
    global w, Gdim, nclus, nmixed
    nw = len(w)
    Hyb = []
    for j in range(nclus):
        d = clusters[j]*nmixed
        Hyb.append(np.zeros((nw, d, d), dtype=np.complex128))

    for i in range(nw):
        for j in range(nclus):
            Hyb[j][i, :, :] = pyqcm.hybridization_function(j, w[i], spin_down)

    return Hyb

