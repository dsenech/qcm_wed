import numpy as np
import matplotlib.pyplot as plt
import re

def draw_operator(file, op_name, show_labels=True, show_neighbors=False, values=False):

    fin = open(file, 'r')

    #-------------------------------------------------------------------------
    # reading positions of sites
    L = ''
    while " sites " not in L:
        L = fin.readline()
        if not L:
            print('file ended without sites')

    L = fin.readline()
    sites = []
    band = []
    cluster = []
    while True:
        L = fin.readline()
        if L == '\n': break
        X = re.split("[(,)\t ]+", L)
        # print('site : ',X)
        sites += [(int(X[4]), int(X[5]), int(X[6]))]
        band += [int(X[3])-1]
        cluster += [int(X[1])-1]

    #-------------------------------------------------------------------------
    # reading basis
    L = ''
    while "phys:" not in L:
        L = fin.readline()
        if not L:
            print('file ended without physical basis')

    basis = []
    for i in range(3):
        L = fin.readline()
        X = re.split("[(,) ]+", L)
        basis += [(float(X[1]), float(X[2]), float(X[3]))]
    B = np.array(basis)
    B = np.linalg.inv(B)
    print('basis = \n', B)        

    S = np.zeros((len(sites), 3))
    for i,s in enumerate(sites):
        S[i,:] = (np.array(s)@B).T
    print('sites:\n', S)        
    print('bands:\n', band)

    #-------------------------------------------------------------------------
    # reading superlattice
    L = ''
    while "superlattice:" not in L:
        L = fin.readline()
        if not L:
            print('file ended without superlattice')

    super = []
    L = fin.readline()
    X = re.split(" ", L)
    dimension = int(X[1])

    for i in range(dimension):
        L = fin.readline()
        X = re.split("[(,) ]+", L)
        super += [(float(X[1]), float(X[2]), float(X[3]))]
    print('superlattice = \n', super)        

    #-------------------------------------------------------------------------
    # reading neighbors
    L = ''
    while "neighbors" not in L:
        L = fin.readline()
        if not L:
            print('file ended without neighbors')

    neighbors = []
    while True:
        L = fin.readline()
        if L == '\n': break
        X = re.split("[(,): ]+", L)
        print('N : ', X)
        neighbors += [(float(X[1]), float(X[2]), float(X[3]))]
    Neigh = np.zeros((len(neighbors), 3))
    for i,s in enumerate(neighbors):
        Neigh[i,:] = (np.array(s)@B).T
    print('neighbors = \n', Neigh)        

    #-------------------------------------------------------------------------
    # reading operator

    L = ''
    while "lattice operators" not in L:
        L = fin.readline()
        if not L:
            raise RuntimeError('file ended without specified lattice operators')

    nc = len(op_name)
    while op_name != L[0:nc]:
        L = fin.readline()
        if not L:
            raise RuntimeError(f'file ended without the operator {op_name}')

    X = re.split(r"[\t(): ]+", L)
    op_type = X[1]
    print(X[0], X[1])
    
    fin.readline()
    elements = []
    while True:
        L = fin.readline()
        if L == '\n' or 'elements' in L: break
        X = re.split("[,\;): ]+", L)
        if len(X) == 6:
            complex = True
            X[3] = X[3][1:]
            v = float(X[3]) + float(X[4])*1j
        else:
            X[-1] = X[-1][0:-1]
            v = float(X[3])
        X[0] = X[0][1:]

        if X[0][-1] == '-': continue
        I = int(X[0][0:-1])
        J = int(X[1][0:-1])
        neighbor = int(X[2])
        if neighbor : continue
        if I > J : continue
        elements += [(I-1,J-1,v)]

    print(elements)

    #-------------------------------------------------------------------------
    # plotting the elements
    for e in elements:
        if np.abs(S[e[0],2]-S[e[1],2]) > 0.001:
            pf = 'r--'
        else:
            pf = 'r-'
        if cluster[e[0]] == cluster[e[1]]:
            alpha = 1.0
        else:
            alpha = 0.5
        if np.linalg.norm(S[e[0],0:2] -  S[e[1],0:2]) < 0.0001 :
            plt.plot([S[e[0],0]], [S[e[0],1]], 'ro', ms = 18, c='w', mec='r', mew=2, alpha = alpha)
        else:
            plt.plot([S[e[0],0], S[e[1],0]], [S[e[0],1], S[e[1],1]], pf, mew=2, alpha = alpha)
            if values:
                plt.text(0.5*(S[e[0],0]+S[e[1],0]), 0.5*(S[e[0],1]+S[e[1],1]), f'${np.round(e[2],5)}$', va='bottom', ha='center', c='r')
    
    #-------------------------------------------------------------------------
    # plotting the sites
    bcol = ['k', 'r', 'b', 'g', 'c', 'm', 'y']
    ncol = len(bcol)
    plt.gca().set_aspect(1)
    plt.gca().axis('off')
    offset = 0.03*np.max(S)

    zmin = np.min(S[:,2])
    zmax = np.max(S[:,2])
    for i in range(S.shape[0]):
        if S[i,2]-0.001 < zmin:
            pass
        else:
            plt.plot(S[i,0], S[i,1], 'ko', ms = 6 + 6*(S[i,2]-zmin)/(zmax-zmin+0.1), mfc='w', c=bcol[band[i]%ncol])
            if show_labels: plt.text(S[i,0], S[i,1]+offset, f'${i+1}$', va='bottom', ha='center', color='b', fontsize=8)
    for i in range(S.shape[0]):
        if S[i,2]-0.001 < zmin:
            plt.plot(S[i,0], S[i,1], 'ko', ms = 6, c=bcol[band[i]%ncol])
            if show_labels: plt.text(S[i,0], S[i,1]-offset, f'${i+1}$', va='top', ha='center', color='b', fontsize=8)

    #-------------------------------------------------------------------------
    # plotting the neighbors

    if show_neighbors:
        for j in range(Neigh.shape[0]):
            for i in range(S.shape[0]):
                plt.plot(S[i,0]+Neigh[j,0], S[i,1]+Neigh[j,1], 'ko', ms = 6, alpha=0.5)

    #-------------------------------------------------------------------------
    plt.title(f'${op_name}$', pad=18)
    plt.show()





    
