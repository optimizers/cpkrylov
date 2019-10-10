cpkrylov v. 1.0
===============

Authors
-------

Daniela di Serafino, University of Campania "Luigi Vanvitelli", Caserta, Italy,
daniela.diserafino@unicampania.it.

Dominique Orban, GERAD and École Polytechnique, Montréal, QC, Canada,
dominique.orban@gerad.ca

Last update
-----------

October 8, 2019

Description
-----------

cpkrylov is a Matlab package implementing Constraint-Preconditioned variants of
Krylov solvers for the solution of the regularized saddle-point linear system

   [ A  B' ] [ x1 ] = [ b1 ]
   [ B  -C ] [ x2 ]   [ b2 ],

where A is n x n, B is m x n, and C is symmetric of size m x m, with m <= n.
The saddle-point matrix is assumed to be nonsingular.

In particular, cpkrylov implements constraint-preconditioned variants of the
CG, Lanczos-CG, MINRES, SYMMLQ, GMRES(l) and DQGMRES methods.

Each solver uses a constraint preconditioner of the form

   [ G  B' ]
   [ B  -C ]

where G is a symmetric approximation to A, and must satisfy the
following condition:

let C = EDE' be a decomposition of C with D nonsingular;
the block diagonal matrix F = blkdiag(G, inv(D)) is positive definite
on the nullspace of [B E], i.e., z'Gz > 0 for all nonzero z in the
nullspace of [B E].

Details on the solvers and the constraint preconditioner are provided in
   Daniela di Serafino and Dominique Orban,
   Constraint-Preconditioned Krylov Solvers for Regularized Saddle-Point Systems,
   Cahier du GERAD G-2019-72, GERAD, Montreal, October 2019.

Note that the argument A may be a matrix or a linear operator, but C, B
and G must be explicit matrices. The preconditioner M is a linear
operator too, such that M*z returns the solution of

   [ G   B' ] [r] = [z1]
   [ B  -C  ] [u]   [z2].

The operator M also implements iterative refinement and residual update,
as suggested (for the case C = 0) in 
   Nicholas I. M. Gould, Mary E. Hribar, and Jorge Nocedal,
   On the solution of equality constrained quadratic programming
   problems quadratic programming problems arising in optimization,
   SIAM Journal on Scientific Computing, 23(4), pp. 1376?1395, 2001.

The linear operatos are defined using the Spot Toolbox by Ewout van
den Berg and Michael P. Friedlander. See
   https://github.com/mpf/spot

Installation
------------

cpkrylov runs under MATLAB (it has been tested under MATLAB 2018b) and requires
the Spot Toolbox (see https://github.com/mpf/spot). In order to use Spot, the
spot-master directory must be added to the Matlab path, e.g., with the Matlab command
   addpath path-to-spot-master

The user must also add to the MATLAB path the root cpkrylov directory and its subdirectories
kernels, ops and utils. This can be done by running cpk_path_setup.m.

Content of cpkrylov package
---------------------------

- reg_cpkrylov.m:         main driver, which performs pre-processing operations,
                          calls the requested solver, performs post-processing operations,
                          and returns solutions and statistics to the user;
- kernels/cpcg.m:         function implementing the constraint-preconditioned CG mathod;
- kernels/cpcglanczos.m:  function implementing the constraint-preconditioned Lanczos version of CG;
- kernels/cpminres.m:     function implementing the constraint-preconditioned MINRES method;
- kernels/cpminres.m:     function implementing the constraint-preconditioned MINRES method;
- kernels/cpsymmlq.m:     function implementing the constraint-preconditioned SYMMLQ method;
- kernels/cpdqgmres.m:    function implementing the constraint-preconditioned DQGMRES method;
- kernels/cpgmres.m:      function implementing the constraint-preconditioned GMRES(l) method;
- ops/LDL2.m:             operator representing the LDL factorization of a symmetric indefinite
                          matrix with optional iterative refinement (needed to apply the constraint
                          preconditioner);
- utils/SymGivens:        function implementing a symmetric Given rotation (by M. A. Saunders and
                          S.-C. Choi), called by cpdqgmres and cpgmres.

Examples of use of cpkrylov
---------------------------

- cpk_exprog1.m:          example program 1, runs cpminres (or cpgc, or cpcglanczos, or
                          cpdqgmres) on the symmetric saddle-point linear system stored in
                          cvxqp1_m_2x2_symm_iter10.mat;
- cpk_exprog2.m:          example program 2, runs cpgmgmres (or cpdqgmres) on the nonsymmetric
                          saddle-point linear system stored in cvxqp2_s_3x3_nonsymm_perm_iter10.mat
- cvxqp1_m_2x2_symm_iter10.mat, cvxqp2_s_3x3_nonsymm_perm_iter10.mat:
                          MATLAB mat-files containing the data for the example programs
