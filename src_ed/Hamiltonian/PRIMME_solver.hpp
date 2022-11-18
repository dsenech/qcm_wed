/*
Interface to use PRIMME eigensolver for Exact Diagonalisation
*/

#ifndef PRIMME_solver_h
#define PRIMME_solver_h

#include <iostream>
#include <vector>
#include "global_parameter.hpp"
#include "primme.h"
#include <Eigen/Sparse>

/**
 Wrapper around dprimme and zprimme routine for type double and complex
 */
template<typename HilbertField>
int call_primme(double* evals, HilbertField* evecs, double* rnorm, primme_params* primme);

/**
 Matrix vector multiplication for Eigen Hamiltonian type
 */
template<typename HilbertField>
void PRIMME_Eigen_matmul(
    void *x, 
    int64_t *ldx,
    void *y,
    int64_t *ldy,
    int *blockSize,
    primme_params *primme,
    int *ierr
) {
     Eigen::Map< Eigen::Matrix<HilbertField,Eigen::Dynamic,1> > xe((HilbertField*) x, primme->n);
     Eigen::Map< Eigen::Matrix<HilbertField,Eigen::Dynamic,1> > ye((HilbertField*) y, primme->n);
     ye = *((Eigen::SparseMatrix<double,Eigen::RowMajor>*) primme->matrix) * xe;
     *ierr = 0;
}

/**
 Matrix vector multiplication for CSR Hamiltonian type
 */
template<typename HilbertField>
void PRIMME_CSR_matmul(
    void *x, 
    int64_t *ldx,
    void *y,
    int64_t *ldy,
    int *blockSize,
    primme_params *primme,
    int *ierr
) {
    //TODO
}

/**
 Compute the ground state of the Hamiltonian using PRIMME eigensolver
 */
template<typename T, typename HilbertField>
void PRIMME_state_solver(
    T* H, //hamiltonian in its different format
    const size_t &dim,
    double &eval, //the returned eigenvalue
    std::vector<HilbertField> &evec, //the returned eigenvector
    const bool verb=false
) {
    /* Solver arrays and parameters */
    double rnorm;   /* Array with the computed eigenpairs residual norms */
    primme_params primme; /* PRIMME configuration struct */
    
    /* Other miscellaneous items */
    int ret;
    int i;

    /* Set default values in PRIMME configuration struct */
    primme_initialize(&primme);
  
   /* Set problem parameters */
   primme.n = (long int) dim; /* set problem dimension */
   primme.numEvals = 1;   /* Number of wanted eigenpairs */
   primme.eps = global_double("accur_lanczos"); /* ||r|| <= eps * ||matrix|| */
   primme.target = primme_smallest;  /* Wanted the smallest eigenvalues */
   primme.maxMatvecs = global_int("max_iter_lanczos");
   if (verb) {
       primme.outputFile = stdout;
       primme.printLevel = 5;
   }
   primme.maxBasisSize = 14;
   primme.minRestartSize = 4;
   primme.maxBlockSize = 1;

    /* Set problem matrix */
    if (Hamiltonian_format == H_format_eigen) {
        primme.matrix = H->H_ptr;
        primme.matrixMatvec = PRIMME_Eigen_matmul<HilbertField>;
    }
    else {
        qcm_ED_throw("Diagonalisation of Hamiltonian of chosen type not implemented with PRIMME eigensolver");
    }
    
   /* Set preconditioner (optional) */
   //primme.applyPreconditioner = LaplacianApplyPreconditioner;
   //primme.correctionParams.precondition = 1;

   /* Set method to solve the problem */
   primme_set_method((primme_preset_method) global_int("PRIMME_algorithm"), &primme);

   /* Call primme */
   ret = call_primme(&eval, evec.data(), &rnorm, &primme);
   
   
   /* Reporting (optional) */
   if (verb) {
       std::cout << "Eval: " << eval << ", rnorm: " << rnorm << std::endl;
       std::cout << "Tolerance : " << primme.aNorm*primme.eps << std::endl;
       std::cout << "Iterations : " << primme.stats.numOuterIterations << std::endl;
       std::cout << "Restarts : " << primme.stats.numRestarts << std::endl;
       std::cout << "Matvecs : " << primme.stats.numMatvecs << std::endl;
       std::cout << "Preconds : " << primme.stats.numPreconds << std::endl;
       if (primme.stats.lockingIssue) {
           std::cout << "A locking problem has occurred" << std::endl;
           std::cout << "Some eigenpairs do not have a residual norm less than the tolerance." << std::endl;
           std::cout << "However, the subspace of evecs is accurate to the required tolerance." << std::endl;
       }
   }

   if (ret != 0) {
      std::cout << "Error: primme returned with nonzero exit status: " << ret << std::endl;
   }
}


#endif
