############
Introduction
############

What is pyqcm?
==============

Pyqcm is a python module that interfaces with two libraries written in C++: **qcm** and **qcm_ED**.
These two libraries provide a collection of functions that help implement quantum cluster methods.
Specifically, qcm_ED provides an exact diagonalization solver for small clusters on which a Hubbard-like model is defined, whereas qcm provides the functions to define infinite-lattice models and to embed the clusters into the lattice via *Cluster Pertrubation Theory* (CPT). Methods like the *Variational Cluster Approximation* (VCA) and *Cluster Dynamical Mean Field Theory* (CDMFT) are then implemented from qcm by the pyqcm module, which is written in Python only.

This document is not a review of the above methods; the reader is referred to the appropriate review articles for that. It is strictly a user's manual for the pyqcm module.
Some degree of understanding of the above methods is however necessary to proceed.
This document also provides insights as to the inner working of the code, and therefore consitutes also a embryonic developer's manual.

Requirements
============

Pyqcm is written in Python and thus requires no compilation.
However, the following librairies are needed:

**qcm**. The source code can be cloned with the following command::

    git clone https://dsenech@bitbucket.org/dsenech/qcm.git

It is written in C++, and requires the qcm_ED library, lapack and python3 for interfacing.

**qcm_ED**. The source code can be cloned with the following command::
    
    git clone https://dsenech@bitbucket.org/dsenech/qcm_ed.git

It is written in C++. It requires lapack and python3 for interfacing.

**CUBA** (http://www.feynarts.de/cuba/) is required by qcm for performing integrals in dimension > 2 (hence the name, for *cubature*). qcm uses version 4.0 or above. Version 2.1 is not supported. However, the CUBA library must be compiled with a modified makefile in order to use it within another shared library (the option -fPIC must be added. See the file `INSTALL` included in the distribution for details).

**BLAS-LAPACK**. Needed for efficient vector and matrix operations.

Compiling qcm and qcm_ED
========================

Please consult the INSTALL files in each library's distribution to see how to compile and install the libraries.
The shared object files qcm.so and qcm_ED.so should be somewhere in the system's Python path.

