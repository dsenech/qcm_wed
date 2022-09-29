import pyqcm
import os
import subprocess
import hopping
#os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np


def dvmc_solver(nproc=1,gs_only=False,exc_only=False, read_soln=False, rand_exc=False, 
                read_gs=True, tol=1e-10):

    nthreads = os.getenv('SLURM_CPUS_PER_TASK')
    nproc    = int(os.getenv('SLURM_NTASKS'))

    print('nthreads = ', nthreads)
    print('nproc = ', nproc)
    
    os.environ["OMP_NUM_THREADS"] = nthreads
    
    print('Printing hopping matrix')
    # printing the model data in a file
    hop = pyqcm.cluster_hopping_matrix(full=True)
    np.save('hop.npy',hop)
    print('Printing interaction')
    I = np.array(pyqcm.interactions())
    np.save('interaction.npy',I)

    """
    for i in range(hop.shape[0]):
        for j in range(hop.shape[1]):
            print('%5.4f '%(hop[i,j].real),end=" ")
        print()
    """

    paramfile = open('params','r')
    for line in paramfile:
        param_name = line.split()[0]
        param_val = line.split()[1]
        if(param_name == 'Nb'):
            Nb = int(param_val)
        if(param_name == 'Nc'):
            Nc = int(param_val)

    Ns = Nb + Nc
    #sec = ' '
    #print("sector: ", sec)
    print('reading sector')
    f = open('sec','r')
    sec = f.readline()
    f.close()
    # print parameters to write solution for QCM
    f = open('QCM_params.def', 'w')
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
    
    if(os.path.exists('model.py')):
        f = open('model.py','r')
    else:
        print('model file not found')
        exit()

    operator_list=[]
    keys = param_set.keys()
    for key in keys:
        if('_1' in key):
            operator_list.append(key[:-2])
    
    operator_list.remove('U')
            
    hopping.write_hop_file(operator_list,'model.py',Ns)

    print("Done writing hop.dat")
    
    p = subprocess.run(["init_params.py", "params"])
    if(p.returncode!=0):
        print("Param files not written properly")
        exit()
    #os.system("init_params.py params")

    # calling the dvmc solver
    if(not exc_only and not read_soln):
        print('Beginning ground state VMC calculation')

        # ground state calculation
        if(read_gs):
            #run_cmd = "mpirun -n "+repr(nproc)+" dvmc.out namelist.def ./output/zqp_opt.dat"
            p = subprocess.run(["mpirun","-n",repr(nproc),"dvmc.out","namelist.def","./output/zqp_opt.dat"])
            #p = subprocess.run(["mpirun","dvmc.out","namelist.def","./output/zqp_opt.dat"])
        else:
            #run_cmd = "mpirun -n "+repr(nproc)+" dvmc.out namelist.def"
            p = subprocess.run(["mpirun","-n",repr(nproc),"dvmc.out","namelist.def"])
            #p = subprocess.run(["mpirun","dvmc.out","namelist.def"])
            #p = subprocess.Popen(["mpirun","dvmc.out","namelist.def"])
            #p.wait()
    if(p.returncode!=0):
        print("VMC GS calculation failed")
        exit()

    E0_VMC, E0_err_VMC = np.loadtxt('./output/zqp_opt.dat',usecols=(0,2))

    pyqcm.cdmft.E0_VMC = E0_VMC
    pyqcm.cdmft.E0_err_VMC = E0_err_VMC
        
    if(not gs_only and not read_soln):
        # make excitation list
        p = subprocess.run(["makeExcitation_from_hopping_only_t.py"])
        if(p.returncode!=0):
            print("Failed to write excitation.def")
            exit()
        
        #if(rand_exc):
        #    os.system("randomize_excitations.py")
        #    os.system("mv excitation_new.def excitation.def")

        print('Beginning dynamical VMC calculation')

        # dynamical calculation
        p = subprocess.run(["mpirun","-n",repr(nproc),"dvmc.out","namelist_G.def","./output/zqp_opt.dat"])
        #p = subprocess.run(["mpirun","dvmc.out","namelist_G.def","./output/zqp_opt.dat"])
        #p = subprocess.Popen(["mpirun","dvmc.out","namelist_G.def", "./output/zqp_opt.dat"])
        #p.wait()
        
        if(p.returncode!=0):
            print("Calculation of excitations failed")
            exit()

        import glob
        files = glob.glob('output/zvo_nCHAm_nAHCm_0*bin')
        # merge binary outputs
        p = subprocess.run(["mergeOutputBin.py"]+files)

        print('Calculating and printing Qmatrix')
        if(p.returncode!=0):
            print("Failed to merge binary files")
            exit()
        
        # compute and output Green's functions
        p = subprocess.run(["dvmc_spectrum_eigh_w_sqrtS.py", "spectrumpara.def", "output",repr(tol)])
        if(p.returncode!=0):
            print("Failed to write Q-matrix")
            exit()

    # reading the solution
    # here, a file qmatrix.def is read. It contains information about the Hilbert space sector, the ground state energy, and the Q matrix 
    if(not gs_only):
        print('Reading dVMC solution')

        with open('output/qmatrix.def') as f:
            solution = f.read()
        f.close()
        pyqcm.read_cluster_model_instance(solution, 0)

    print('Done reading dVMC solution')
    #subprocess.run(["mv","hop.npy","hop_prev_iter.npy"])
