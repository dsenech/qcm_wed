import pyqcm
import os
#os.environ["OMP_NUM_THREADS"] = "8"
import numpy as np


def dvmc_solver(nproc=1,gs_only=False,exc_only=False, read_soln=False, rand_exc=False, 
                read_gs=True, tol=1e-5):

    # printing the model data in a file
    hop = pyqcm.cluster_hopping_matrix(full=True)
    np.save('hop.npy',hop)
    I = np.array(pyqcm.interactions())
    np.save('interaction.npy',I)

    a,b=pyqcm.properties()
    a = a.split()
    b = b.split()
    sec_ind = a.index('sector_1')
    sec = b[sec_ind]
    print("sector: ", sec)
    # print parameters to write solution for QCM
    f = open('QCM_params.def', 'w')
    f.write('sector '+sec+'\n')
    param_set = pyqcm.parameter_set()
    f.write('sector '+sec+'\n')
    for name, val in param_set.items():
        # print non-zero parameters
        if(abs(val[0])>1e-8):
            if('b' in name):
                # print bath params without underscore
                f.write(name.partition('_')[0]+'  '+repr(val[0])+'\n')
            elif('_' not in name):
                # print cluster params
                f.write(name+'  '+repr(val[0])+'\n')
    f.close()

    os.system("init_params.py params")

    # calling the dvmc solver
    if(not exc_only and not read_soln):
        print('Beginning ground state VMC calculation')

        # ground state calculation
        if(read_gs):
            run_cmd = "mpirun -n "+repr(nproc)+" dvmc.out namelist.def ./output/zqp_opt.dat"
        else:
            run_cmd = "mpirun -n "+repr(nproc)+" dvmc.out namelist.def"
        #os.system("mpirun -n dvmc.out namelist.def")
    
        os.system(run_cmd)


    if(not gs_only and not read_soln):
        # make excitation list
        os.system("makeExcitation_from_hopping_only_t.py")

        if(rand_exc):
            os.system("randomize_excitations.py")
            os.system("mv excitation_new.def excitation.def")

        print('Beginning dynamical VMC calculation')

        # dynamical calculation
        #os.system("dvmc.out namelist_G.def ./output/zqp_opt.dat")
        run_cmd = "mpirun -n "+repr(nproc)+" dvmc.out namelist_G.def ./output/zqp_opt.dat"
        os.system(run_cmd)

        # merge binary outputs
        os.system("mergeOutputBin.py output/zvo_nCHAm_nAHCm_0*bin")

        print('Calculating and printing Qmatrix')

        # compute and output Green's functions
        os.system("dvmc_spectrum_eigh_w_sqrtS.py spectrumpara.def output "+repr(tol))

    # reading the solution
    # here, a file qmatrix.def is read. It contains information about the Hilbert space sector, the ground state energy, and the Q matrix 
    if(not gs_only):
        print('Reading dVMC solution')

        with open('output/qmatrix.def') as f:
            solution = f.read()
        f.close()
        pyqcm.read_cluster_model_instance(solution, 0)
