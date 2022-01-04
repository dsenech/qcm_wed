import numpy as np
import pyqcm
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import hsv_to_rgb


def __color(z):
    x = np.modf(np.angle(z) / (2 * np.pi) + 1.66667)[0]
    return hsv_to_rgb(np.array([x, 0.8, 0.8]))


def plot_profile(n_scale=1, bond_scale=1, current_scale=1, spin_scale=1,
                 spin_angle=0, bond_spin_scale=1, singlet_scale=1, triplet_scale=1, file=None, layer=0):
    """Produces a figure of various local quantities on the repeated unit, from averages computed in the ground state wavefunction

    :param float n_scale: scale factor applied to the density
    :param float bond_scale: scale factor applied to the bond charge density
    :param float current_scale: scale factor applied to the currents on the bonds
    :param float spin_scale: scale factor applied to the spins on the sites
    :param float spin_angle: angle at which to draw the spins, from their nominal direction
    :param float bond_spin_scale: scale factor applied to the spins on the bonds 
    :param float singlet_scale: scale factor applied to the singlet pairing amplitudes
    :param float triplet_scale: scale factor applied to the triplet pairing amplitudes
    :param str file: name of the output file, if not None
    :param float z: layer number (z coordinate)
        
    """


    S, B = pyqcm.site_and_bond_profile()

    # printing on the screen

    print('\nsite profiles:')
    for x in S:
        # if np.abs(x[2] - z) > 0.01 :
        #     continue
        print('({: 1.1f},{: 1.1f},{: 1.1f}) : n={: 1.4f}    S=({: 1.4f},{: 1.4f})    psi={: 1.4f}'.format(
            x[0], x[1], x[2], x[3], x[4], x[6], x[7] + x[8] * 1j))

    print('\nbond profiles:')
    for x in B:
        # if np.abs(x[2] - z) > 0.01 :
        #     continue
        print("""({: 0.1f},{: 0.1f},{: 1.1f})-({: 0.1f},{: 0.1f},{: 0.1f}) : b0={: 1.3f}    b=({: 1.3f},{: 1.3f},{: 1.3f})\n
                d0={: 1.3f}   d=({: 1.3f},{: 1.3f},{: 1.3f})""".format(
            x[0].real, x[1].real, x[2].real, x[0].imag, x[1].imag, x[2].imag, x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10]))

    # plotting

    plt.subplot(221)
    ax = plt.gca()
    ax.set_aspect(1)
    plt.xlim(np.ma.min(S[:, 0]) - 0.5, np.ma.max(S[:, 0]) + 0.5)
    plt.ylim(np.ma.min(S[:, 1]) - 0.5, np.ma.max(S[:, 1]) + 0.5)

    ns = len(S)
    nb = len(B)
    bond_scale *= 0.7
    t = spin_angle*np.pi/180
    rotation_matrix = np.array([[np.cos(t), np.sin(t)], [-np.sin(t), np.cos(t)]])

    #--------------------------------------------------------------
    # normal part

    plt.title('charge and spin', fontsize=9)
    plt.plot(S[:,0], S[:,1], 'ko', ms = 1)

    for i in range(nb):
        if np.abs(B[i,2].real - layer) > 0.01 or np.abs(B[i,2].imag - layer):
            continue
        plt.plot([B[i, 0].real, B[i, 0].imag], [B[i, 1].real, B[i, 1].imag], 'k-', lw=0.5)
        x = 0.5 * (B[i, 0].real + B[i, 0].imag)
        y = 0.5 * (B[i, 1].real + B[i, 1].imag)
        dx = B[i, 0].real - B[i, 0].imag
        dy = B[i, 1].real - B[i, 1].imag
        v = B[i, 3].real * bond_scale
        current = B[i, 3].imag * current_scale
        spin = np.array([B[i, 4].real, B[i, 6].real])
        spin = np.dot(rotation_matrix, spin) * bond_spin_scale
        sx = spin[0]
        sz = spin[1]

        if v > 0:
            ax.add_patch(patches.Ellipse((x,y), v, 0.3*v, angle = np.degrees(np.arctan2(dy,dx)), facecolor='blue', alpha=0.4))
        else:
            ax.add_patch(patches.Ellipse((x,y), v, 0.3*v, angle = np.degrees(np.arctan2(dy,dx)), facecolor='red', alpha=0.4))
        if np.abs(current) > 0.01 :
            ax.add_patch(patches.Arrow(x-0.5*current*dx, y-0.5*current*dy, current*dx, current*dy, width=0.2, fc = 'red', ec = 'black', lw = 0.5))
        if np.sqrt(sx*sx+sz*sz)>0.01 :
            ax.add_patch(patches.Arrow(x-0.5*sx, y-0.5*sz, sx, sz, width=0.2, fc = 'orange', ec = 'black', lw = 0.5))

    for i in range(ns):
        if np.abs(S[i,2] - layer) > 0.01 :
            continue
        x = S[i, 0]
        y = S[i, 1]
        n = S[i, 3]-1
        spin = np.array([S[i, 4], S[i, 6]])
        spin = np.dot(rotation_matrix,spin)*spin_scale
        sx = spin[0]
        sz = spin[1]
        if n>0 :
            ax.add_patch(patches.Circle((x,y), n_scale*n, facecolor='blue', alpha=0.4))
        else:
            ax.add_patch(patches.Circle((x,y), -n_scale*n, facecolor='red', alpha=0.4))
        if np.sqrt(sx*sx+sz*sz)>0.01 :
            ax.add_patch(patches.Arrow(x-0.5*sx, y-0.5*sz, sx, sz, width=0.1, fc = 'red', ec = 'black', lw = 0.5, zorder=100))

    #--------------------------------------------------------------
    # singlet part

    plt.subplot(222)
    ax = plt.gca()
    ax.set_aspect(1)
    plt.xlim(np.ma.min(S[:,0]) - 1, np.ma.max(S[:,0]) + 1)
    plt.ylim(np.ma.min(S[:,1]) - 1, np.ma.max(S[:,1]) + 1)

    plt.title('singlet pairing', fontsize=9)
    plt.plot(S[:,0], S[:,1], 'ko', ms = 1)

    for i in range(nb):
        if np.abs(B[i,2].real - layer) > 0.01 or np.abs(B[i,2].imag - layer):
            continue
        plt.plot([B[i, 0].real, B[i, 0].imag], [B[i, 1].real, B[i, 1].imag], 'k-', lw = 0.5)
        x = 0.5*(B[i, 0].real+B[i, 0].imag)
        y = 0.5*(B[i, 1].real+B[i, 1].imag)
        dx = B[i, 0].real-B[i, 0].imag
        dy = B[i, 1].real-B[i, 1].imag
        v = np.abs(B[i, 7])*singlet_scale
        col = __color(B[i, 7])
        ax.add_patch(patches.Ellipse((x,y), v, 0.3*v, angle = np.degrees(np.arctan2(dy,dx)), facecolor=col, alpha=0.4))

    for i in range(ns):
        if np.abs(S[i,2] - layer) > 0.01 :
            continue
        x = S[i, 0]
        y = S[i, 1]
        psi = S[i, 7]+S[i, 8]*1.0j
        v = np.abs(psi)*singlet_scale
        col = __color(psi)
        ax.add_patch(patches.Circle((x,y), v, facecolor=col, alpha=0.4))

    #--------------------------------------------------------------
    # dz part

    plt.subplot(223)
    ax = plt.gca()
    ax.set_aspect(1)
    plt.xlim(np.ma.min(S[:,0]) - 1, np.ma.max(S[:,0]) + 1)
    plt.ylim(np.ma.min(S[:,1]) - 1, np.ma.max(S[:,1]) + 1)

    plt.title('triplet pairing (d_z)', fontsize=9)
    plt.plot(S[:,0], S[:,1], 'ko', ms = 1)

    for i in range(nb):
        if np.abs(B[i,2].real - layer) > 0.01 or np.abs(B[i,2].imag - layer):
            continue
        plt.plot([B[i, 0].real, B[i, 0].imag], [B[i, 1].real, B[i, 1].imag], 'k-', lw = 0.5)
        x = 0.5*(B[i, 0].real+B[i, 0].imag)
        y = 0.5*(B[i, 1].real+B[i, 1].imag)
        dx = B[i, 0].real-B[i, 0].imag
        dy = B[i, 1].real-B[i, 1].imag
        v = np.abs(B[i, 10])*triplet_scale
        z = B[i, 10]
        if np.abs(np.angle(z)) > np.pi/2:
            z = -z
            v *= -1
        col = __color(z)
#       ax.add_patch(patches.Ellipse((x,y), v, 0.3*v, angle = np.degrees(np.arctan2(dy,dx)), facecolor=col, alpha=0.4))
        ax.add_patch(patches.Arrow(x-0.5*v*dx, y-0.5*v*dy, v*dx, v*dy, width=0.5, fc = col, ec = 'black', lw = 0.5, zorder=100))

    #--------------------------------------------------------------
    # dx, dy part

    plt.subplot(224)
    ax = plt.gca()
    ax.set_aspect(1)
    plt.xlim(np.ma.min(S[:,0]) - 1, np.ma.max(S[:,0]) + 1)
    plt.ylim(np.ma.min(S[:,1]) - 1, np.ma.max(S[:,1]) + 1)

    plt.title('triplet pairing (d_x and d_y)', fontsize=9)
    plt.plot(S[:,0], S[:,1], 'ko', ms = 1)

    for i in range(nb):
        if np.abs(B[i,2].real - layer) > 0.01 or np.abs(B[i,2].imag - layer):
            continue
        plt.plot([B[i, 0].real, B[i, 0].imag], [B[i, 1].real, B[i, 1].imag], 'k-', lw = 0.5)
        x = 0.5*(B[i, 0].real+B[i, 0].imag)
        y = 0.5*(B[i, 1].real+B[i, 1].imag)
        dx = B[i, 0].real-B[i, 0].imag
        dy = B[i, 1].real-B[i, 1].imag
        # spin = np.array([np.abs(B[i, 8]), np.abs(B[i, 9])])
        # spin = np.dot(rotation_matrix,spin)*triplet_scale
        # sx = spin[0]
        # sz = spin[1]
        # col = __color(B[i, 8]/(B[i, 9]+1e-6))
        # if np.sqrt(sx*sx+sz*sz) > 0.01 :
        #   ax.add_patch(patches.Arrow(x-0.5*sx, y-0.5*sz, sx, sz, width=0.2, fc = col, ec = 'black', lw = 0.5, zorder=100))
        spinR = np.array([B[i, 8].real, B[i, 9].real])
        spinI = np.array([B[i, 8].imag, B[i, 9].imag])
        spinR = np.dot(rotation_matrix,spinR)*triplet_scale
        spinI = np.dot(rotation_matrix,spinI)*triplet_scale
        sx = spinR[0]
        sz = spinR[1]
        if np.sqrt(sx*sx+sz*sz) > 0.01 :
            ax.add_patch(patches.Arrow(x-0.5*sx, y-0.5*sz, sx, sz, width=0.2, fc = 'blue', ec = 'black', lw = 0.5, zorder=100))
        sx = spinI[0]
        sz = spinI[1]
        if np.sqrt(sx*sx+sz*sz) > 0.01 :
            ax.add_patch(patches.Arrow(x-0.5*sx, y-0.5*sz, sx, sz, width=0.2, fc = 'red', ec = 'black', lw = 0.5, zorder=100))

    plt.gcf().suptitle(pyqcm.parameter_string()+', layer '+str(layer), y=1.0, fontsize=9)
    plt.tight_layout()

    if file is None:
        plt.show()
    else:
        plt.savefig(file)
