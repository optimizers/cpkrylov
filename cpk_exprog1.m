% CPKRYLOV EXAMPLE PROGRAM 1
%
%==========================================================================
% August 23, 2019.
% Daniela di Serafino, daniela.diserafino@unicampania.it.
% Dominique Orban, dominique.orban@gerad.ca.
%
%==========================================================================
%
% This program runs CP-MINRES (CP-CG, CP-CGLANCZOS, CP-DQGMRES) on a
% saddle-point linear system whose matrix has the form
%
%                     [ H+rho*I        0        |    J1'   ]     
%       [ Q |  A']    [    0     S^(-1)*Z+rho*I |    J2'   ]
%   K = ----------  = --------------------------|----------- ,
%       [ A | -C ]    [   J1           J2       | -delta*I ]
%
% where
% - Q is symmetric and n x n,
% - A is m x n, with m <= n,
% - C is m x m.
%
% The matrix satisfies additional conditions, as explained in
%   D. di Serafino and D. Orban,
%   Constraint-Preconditioned Krylov Solvers for Regularized
%   Saddle-Point Systems.
%   TBA
%
% The linear system comes from a collection of regularized saddle-point
% systems, generated by applying a primal-dual regularized interior point
% method to quadratic programming problems from CUTEst.
% See
%   D. Orban,
%   A collection of linear systems arising from interior-point methods for
%   quadratic optimization,
%   Cahier du GERAD G-2015-117, GERAD, Montr�al, QC, Canada,
% for a description of the collection and details on the matrix K reported
% above.
%
%==========================================================================

% Load saddle-point system (problem cvxqp1-m, interior point iter 10)
filename = 'cvxqp1_m_2x2_symm_iter10.mat';
load(filename);               % load K, n=dim(K), nH, nJ, nZ, rhs
n = nH;
m = nJ;
clear nH nJ nZ;

fprintf('\n\n==================================================================');
fprintf('\n                   cpkrylov example program 1');
fprintf('\n==================================================================\n');
fprintf('\nsymmetric saddle-point system from');
fprintf('\n- quadratic programming problem cvxqp1-m (n = %d, m = %d)', n, m);
fprintf('\n- interior point iteration 10\n');

% Build constraint preconditioner
Q = K(1:n,1:n);
dd = spdiags(Q, 0);
G = spdiags(dd, 0, n, n);
A = K(n+1:end,1:n);
C = - K(n+1:end,n+1:end);
CP = [G  A'; A -C];
    
% Choose solver
cpk = @cpminres;
cpkstring = 'CP-MINRES';
% cpk = @cpcg;
% cpkstring = 'CP-CG';
% cpk = @cpcglanczos;
% cpkstring = 'CP-CGLANCZOS';
% cpk = @cpdqgmres;  opts.mem = 2;
% cpkstring = 'CP-DQGMRES(2)';
fprintf('\n**************************** %s ***************************\n', cpkstring);
        
% Set solver options
% Some options are set to default values -- we can avoid their explicit setup
opts.print = false;     % display info during execution;
opts.atol = 1.0e-6;     % abs tol for stop criterion (default)
opts.rtol = 1.0e-6;     % rel tol for stop criterion (default)
opts.itmax = 500;      % max num solver iterations
fprintf('\natol = %8.2e,  rtol = %8.2e,  itmax = %d', ...
    opts.atol, opts.rtol, opts.itmax);

% Set options for residual update and iterative refinement
opts.residual_update = true; % use residual update (true/false)
opts.nitref = 1;             % max # iterative refinement steps
opts.force_itref = true;     % force iterative refinement (true/false)
opts.itref_tol = 1.0e-8;     % rel tol for iterative refinement (default)
fprintf('\nresidual_update = %d,  nitref = %d,  force_iref = %d,  itref_tol = %7.1e\n', ...
    opts.residual_update, opts.nitref, opts.force_itref, opts.itref_tol);

    
% Call CP-Krylov solver
ts = tic;
[cpkx,  cpkstats, cpkflag] = reg_cpkrylov(cpk, rhs, Q, A, C, G, opts);
ttot = toc(ts);

% Solve K*x = rhs with a direct method (for comparison)
x = K \ rhs;

% Print error and other info
fprintf('\n2-norm relative error in the solution = %8.2e', norm(x-cpkx)/norm(x));
fprintf('\niters = %d,  solved (1 yes, 0 no) = %d', cpkstats.niters, cpkflag.solved);
fprintf('\ntime (prec setup, solve, reg_cpkrylov) = %9.3e,  %9.3e,  %9.3e\n\n', ...
        cpkstats.ptime, cpkstats.stime, ttot);
    
% Print residual norm history
fig = figure;
semilogy(0: size(cpkstats.residHistory,1)-1, cpkstats.residHistory, 'b-.', ...
    'LineWidth', 2);
legend(cpkstring, 'fontsize', 14);
xlabel('iters');
ylabel('||res||');
title('residual history (cvxqp1-m, ip iter 10, symm matrix)');
set(gca, 'fontsize', 14);

