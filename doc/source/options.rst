.. include:: options_intro.txt

Boolean options
=========================
.. csv-table::
    :header: "name", "default", "description"
    :widths: 15, 10, 50

    "modified_Lanczos", "false", "Uses the modified Lanczos method for the ground state instead of the usual Lanczos method"
    "CSR_sym_store", "false", "stores CSR matrices fully for openMP application"
    "check_lanczos_residual", "false", "checks the Lanczos residual at the end of the eigenvector computation"
    "one_body_solution", "false", "Only solves the one-body part of the problem, for the Green function"
    "zero_dim", "false", "sets the spatial dimension to zero, on any model"
    "SEF_calc", "false", "Calculates the Self-energy functional, even when not performing VCA"
    "print_Hamiltonian", "false", "Prints the Hamiltonian on the screen, if small enough"
    "nosym", "false", "does not take cluster symmetries into account"
    "no_degenerate_BL", "false", "forbids band lanczos to proceed when the eigenstates have degenerate energies"
    "BFGS", "false", "Uses the BFGS formula in the quasi-Newton method, instead of the symmetric rank-1 formula"
    "potential_energy", "false", "Calculates the potential energy with Tr(Sigma*G) when computing averages"
    "print_all", "false", "prints dependent parameters as well"
    "epsilon_algorithm", "false", "uses the epsilon algorithm for prediction instead of polynomials"
    "NR_box", "false", "In the Newton-Raphson algorithm, makes sure the first move is within the steps"
    "periodized_averages", "false", "computes lattice averages using the periodized Green function"
    "periodic", "false", "considers the cluster(s) as periodic"
    "continued_fraction", "false", "Uses the continued fraction solver for the Green function instead of the band Lanczos method"
    "strip_anomalous_self", "false", "sets to zero the anomalous part of the self-energy"
    "epsilon_predictor", "false", "use epsilon algorithm to periodically predict bath parameters"
    "dual_basis", "false", "uses the dual basis for wavevector computations"



Integer-valued options
=========================
.. csv-table::
    :header: "name", "default", "description"
    :widths: 15, 10, 50

    "verbose", "0", "level of verbosity"
    "max_iter_lanczos", "600", "Maximum number of iterations in the Lanczos procedure"
    "max_iter_CF", "400", "Maximum number of iterations in the continuous fraction Lanczos procedure"
    "print_precision", "8", "precision of printed output"
    "GK_min_regions", "8", "minimum number of regions in the Gauss-Kronrod method"
    "seed", "0", "seed of the random number generator"
    "max_iter_QN", "60", "maximum number of iterations in the quasi-Newton method"
    "dim_max_print", "64", "Maximum dimension for printing vectors and matrices"
    "max_dim_full", "256", "Maximum dimension for using full diagonalization"
    "minimize_max_iter", "50000", "maximum number of iterations in the main minimization loop"
    "cuba3D_mineval", "16000", "minimum number of integrand evaluations in CUBA (3D)"
    "max_iter_NR", "30", "maximum number of iterations in the Newton-Raphson method"
    "cuba2D_mineval", "1024", "minimum number of integrand evaluations in CUBA (2D)"
    "minimize_max_calls", "100000", "maximum number of calls to the function being minimized"
    "sim_ann_updates", "10000", "number of updates in the simulated annealing method"
    "max_iter_BL", "600", "Maximum number of iterations in the band Lanczos procedure"
    "Davidson_states", "1", "Number of states requested in the Davidson-Liu algorithm"
    "kgrid_side", "32", "number of wavevectors on the side in a fixed wavevector grid"
    "max_iter_brent", "10000", "maximum number of iterations in the Brent minimization method"



Real-valued options
=========================
.. csv-table::
    :header: "name", "default", "description"
    :widths: 15, 10, 50

    "minimum_weight", "0.01", "minimum weight in the density matrix"
    "sim_ann_damping", "0.99", "damping of temperature in simulated annealing"
    "band_lanczos_minimum_gap", "1e-05", "gap between the lowest two states in BL below which the method fails"
    "accur_lanczos", "1e-12", "tolerance of the Ritz residual estimate in the Lanczos method"
    "accur_Davidson", "1e-05", "maximum norm of residuals in the Davidson-Liu algorithm"
    "accur_SEF", "5e-08", "Accuracy of the Potthoff functional"
    "large_scale", "20", "high-frequency region for imaginary frequency axis integrals"
    "asy_scale", "1", "scale factor for operators in asy plots"
    "SEF_step", "0.0001", "step used in calculating the SEF gradient (with option -gm only)"
    "accur_deflation", "1e-07", "norm below which a vector is deflated in the Band Lanczos method"
    "Hooke_Jeeves_damping", "0.9", "damping of iterations in Hooke-Jeeves minimization"
    "newton_damping", "1", "damping factor in Newton's method"
    "accur_band_lanczos", "1e-12", "energy difference tolerance for stopping the BL process"
    "temperature", "0", "Temperature of the system."
    "accur_continued_fraction", "0.01", "value of beta below which the simple Lanczod process stops"
    "max_newton_gradient", "1", "maximum value of the gradient in Newton's method before it fails"
    "max_newton_diff", "1", "maximum step in Newton's method"
    "accur_OP", "0.0001", "accuracy of lattice averages"
    "cutoff_scale", "1e+12", "high-frequency cutoff in integrals"
    "accur_Q_matrix", "1e-05", "tolerance in the normalization of the Q matrix"
    "Qmatrix_tolerance", "0.0001", "minimum value of a Qmatrix coefficient"
    "small_scale", "0.5", "low-frequency region for imaginary frequency axis integrals"
    "asy_spin_scale", "0.1", "scale factor for spins in asy plots"
    "eta", "0.0001", "value of the imaginary part of the frequency in Chern number/Berry phase computations"



Char-valued options
=========================
.. csv-table::
    :header: "name", "default", "description"
    :widths: 15, 10, 50

    "Hamiltonian_format", "S", "Desired Hamiltonian format: S (CSR matrix), O (individual operators), F (factorized), N (none = on the fly)"
    "periodization", "G", "periodization scheme: G, S, M, C or N (None)"



