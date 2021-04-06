#include <cstring>
#include <fstream>

#include "global_parameter.hpp"

namespace QCM {
  bool is_global_parameter_initialized = false;
  void global_parameter_init()
  {
    if(is_global_parameter_initialized) return;
    is_global_parameter_initialized = true;

    new_global_bool(false, "BFGS","Uses the BFGS formula in the quasi-Newton method, instead of the symmetric rank-1 formula");
    new_global_bool(false, "dual_basis", "uses the dual basis for wavevector computations");
    new_global_bool(false,"epsilon_algorithm","uses the epsilon algorithm for prediction instead of polynomials");
    new_global_bool(false,"epsilon_predictor","use epsilon algorithm to periodically predict bath parameters");
    new_global_bool(false,"NR_box","In the Newton-Raphson algorithm, makes sure the first move is within the steps");
    new_global_bool(false,"periodic","considers the cluster(s) as periodic");
    new_global_bool(false,"periodized_averages","computes lattice averages using the periodized Green function");
    new_global_bool(false,"potential_energy","Calculates the potential energy with Tr(Sigma*G) when computing averages");
    new_global_bool(false,"SEF_calc","Calculates the Self-energy functional, even when not performing VCA");
    new_global_bool(false,"zero_dim","sets the spatial dimension to zero, on any model");
    new_global_bool(false,"print_all","prints dependent parameters as well");

    new_global_double(0.0001, "eta", "value of the imaginary part of the frequency in Chern number/Berry phase computations");
    new_global_double(0.1,"asy_spin_scale","scale factor for spins in asy plots");
    new_global_double(0.5,"Hooke_Jeeves_damping", "damping in the Hooke-Jeeves minimization method");
    new_global_double(0.5,"small_scale", "low-frequency region for imaginary frequency axis integrals");
    new_global_double(0.98,"sim_ann_damping", "damping in the simulated annealing method");
    new_global_double(1,"asy_scale","scale factor for operators in asy plots");
    new_global_double(1.0,"max_newton_diff", "maximum step in Newton's method");
    new_global_double(1.0,"max_newton_gradient", "maximum value of the gradient in Newton's method before it fails");
    new_global_double(1.0,"newton_damping", "damping factor in Newton's method");
    new_global_double(1.0e12, "cutoff_scale", "high-frequency cutoff in integrals");
    new_global_double(1e-4, "accur_OP", "accuracy of lattice averages");
    new_global_double(1e-4,"SEF_step", "step used in calculating the SEF gradient (with option -gm only)");
    new_global_double(20,"large_scale", "high-frequency region for imaginary frequency axis integrals");
    new_global_double(5e-8,"accur_SEF", "Accuracy of the Potthoff functional");
    new_global_double(0.9, "Hooke_Jeeves_damping", "damping of iterations in Hooke-Jeeves minimization");
    new_global_double(0.99, "sim_ann_damping", "damping of temperature in simulated annealing");
    
    new_global_int(0,"seed","seed of the random number generator");
    new_global_int(1000,"sim_ann_updates", "number of updates at fixed T in the simulated annealing method");
    new_global_int(10000,"max_iter_brent","maximum number of iterations in the Brent minimization method");
    new_global_int(100000,"minimize_max_calls", "maximum number of calls to the function being minimized");
    new_global_int(1024,"cuba2D_mineval","minimum number of integrand evaluations in CUBA (2D)");
    new_global_int(16000,"cuba3D_mineval","minimum number of integrand evaluations in CUBA (3D)");
    new_global_int(30,"max_iter_NR","maximum number of iterations in the Newton-Raphson method");
    new_global_int(32,"kgrid_side","number of wavevectors on the side in a fixed wavevector grid");
    new_global_int(50000,"minimize_max_iter", "maximum number of iterations in the main minimization loop");
    new_global_int(60,"max_iter_QN","maximum number of iterations in the quasi-Newton method");
    new_global_int(64,"dim_max_print","Maximum dimension for printing vectors and matrices");
    new_global_int(8, "GK_min_regions","minimum number of regions in the Gauss-Kronrod method");
    new_global_int(8,"print_precision","precision of printed output");
    new_global_int(10000,"sim_ann_updates","number of updates in the simulated annealing method");

    new_global_char('G', "periodization", "periodization scheme: G, S, M, C or N (None)");

    // ED global parameters

    new_global_bool(false,"check_lanczos_residual","checks the Lanczos residual at the end of the eigenvector computation");
    new_global_bool(false,"no_degenerate_BL","forbids band lanczos to proceed when the eigenstates have degenerate energies");
    new_global_bool(false,"nosym", "does not take cluster symmetries into account");
    new_global_bool(false,"one_body_solution","Only solves the one-body part of the problem, for the Green function");
    new_global_bool(false,"print_Hamiltonian","Prints the Hamiltonian on the screen, if small enough");
    new_global_bool(false,"CSR_sym_store","stores CSR matrices fully for openMP application");
    new_global_bool(false,"strip_anomalous_self","sets to zero the anomalous part of the self-energy");
    new_global_bool(false,"modified_Lanczos","Uses the modified Lanczos method for the ground state instead of the usual Lanczos method");
    new_global_bool(false,"continued_fraction","Uses the continued fraction solver for the Green function instead of the band Lanczos method");


    new_global_double(1e-12,"accur_band_lanczos","energy difference tolerance for stopping the BL process");
    new_global_double(0.01,"accur_continued_fraction","value of beta below which the simple Lanczod process stops");
    new_global_double(1.0e-5,"accur_Davidson","maximum norm of residuals in the Davidson-Liu algorithm");
    new_global_double(1e-7,"accur_deflation","norm below which a vector is deflated in the Band Lanczos method");
    new_global_double(1e-12,"accur_lanczos","tolerance of the Ritz residual estimate in the Lanczos method");
    new_global_double(1.0e-5,"accur_Q_matrix","tolerance in the normalization of the Q matrix");
    new_global_double(1e-5,"band_lanczos_minimum_gap","gap between the lowest two states in BL below which the method fails");
    new_global_double(0.01,"minimum_weight","minimum weight in the density matrix");
    new_global_double(1.0e-4, "Qmatrix_tolerance", "minimum value of a Qmatrix coefficient");
    new_global_double(0.0,"temperature", "Temperature of the system.");

    new_global_int(1,"Davidson_states","Number of states requested in the Davidson-Liu algorithm");
    new_global_int(64,"dim_max_print","Maximum dimension for printing vectors and matrices");
    new_global_int(256,"max_dim_full","Maximum dimension for using full diagonalization");
    new_global_int(600,"max_iter_BL","Maximum number of iterations in the band Lanczos procedure");
    new_global_int(400,"max_iter_CF","Maximum number of iterations in the continuous fraction Lanczos procedure");
    new_global_int(600,"max_iter_lanczos","Maximum number of iterations in the Lanczos procedure");
    new_global_int(0,"seed","seed of the random number generator");
    new_global_int(0,"verbose","level of verbosity");

    new_global_char('S', "Hamiltonian_format", "Desired Hamiltonian format: S (CSR matrix), O (individual operators), F (factorized), N (none = on the fly)");

  }
}

