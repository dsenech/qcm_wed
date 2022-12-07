/*
 Implementation of the Hamiltonian using PETSc for multinode ED
*/

#ifndef Hamiltonian_petsc
#define Hamiltonian_petsc

#include "Hamiltonian_base.hpp"
#include <petscmat.h>

class Hamiltonian_PETSc: public Hamiltonian
{
    Mat H_petsc;
    
    Hamiltonian() {
        PetscInitialize();
        MatCreate(PETSC_COMM_WORLD, &csr);
        MatSetSizes(csr,PETSC_DECIDE,PETSC_DECIDE,n,n);
        MatSetType(csr,MATAIJ);
    }
    
    void populate_Hamiltonian();
    void update_Hamiltonian();
};





Hamiltonian::populate_Hamiltonian() {
    MatSetValues(H_petsc, , ADD_VALUES);
    MatAssemblyBegin(H_petsc, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(H_petsc, MAT_FINAL_ASSEMBLY);
}

Hamiltonian::update_Hamiltonian() {
    MatSetValues(H_petsc, , INSERT_VALUES);
    MatAssemblyBegin(H_petsc, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(H_petsc, MAT_FINAL_ASSEMBLY);
}

#endif
