function [x, y, flags, stats] = cpcg(b, A, C, M, opts)

%======================================================================
% [x, y, flags, stats] = cpcg(b, A, C, M, opts)
%
% Constraint-preconditioned CG for generalized saddle-point systems.
%
%======================================================================
% Last update, September 9, 2017.
% Daniela di Serafino, daniela.diserafino@unicampania.it.
% Dominique Orban, dominique.orban@gerad.ca.
%
%======================================================================
% This function solves the regularized saddle-point system
%
%  [ A   B' ] [x] = [b]
%  [ B  -C  ] [y]   [0],
%
% where A is n x n, B is m x n, C is m x m, with m <= n, and A and C
% are symmetric. The system must also satisfy the following condition:
%
% let C = EDE' be a decomposition of C with D nonsingular;
% the block-diagonal matrix blkdiag(A, inv(D)) must be positive
% definite on the nullspace of [B E].
%
% The method uses a constraint preconditioner of the form
%
%  [ G   B' ]
%  [ B  -C  ],
%
% where G is a symmetric approximation to H, and must be chosen so that
% blkdiag(G, inv(D)) is positive definite on the nullspace of [B E].
%
% The iterations stop when
%
%   (residNorm <= stopTol = atol + rtol * residNorm0)  or  (itn = itmax),
% 
% where residNorm and residNorm0 are the 2-norms of the current and
% initial residuals, atol and rtol are absolute and relative tolerances,
% itn is the iteration index, and itmax is the maximum number of
% iterations.
%
% NOTE that
% - the argument A may be a matrix or a linear operator, but C and G
%   must be explicit matrices;
% - B is not explicitly passed to cpgmres as an argument, but it
%   has been used to form the constraint preconditioner stored
%   in M (see reg_cpkrylov.m);
% - M must be an operator such that M*z returns the solution of
%
%  [ G   B' ] [r] = [z1]
%  [ B  -C  ] [u]   [z2].
%
% The linear operatos are defined using the Spot Toolbox by Ewout van
% den Berg and Michael P. Friedlander.
% See http://www.cs.ubc.ca/labs/scl/spot.
% 
%======================================================================
% REFERENCE
%   H.S. Dollar, N.I.M. Gould, W.H.A. Schilders, and A.J. Wathen,
%   Implicit-Factorization Preconditioning and Iterative Solvers for
%   Regularized Saddle-Point Systems, SIAM Journal on Matrix Analysis
%   and Applications, 28(1), pp. 170-189, 2006.
%
%======================================================================
% INPUT ARGUMENTS:
% b:     n-vector, the vector b in the rhs of the saddle-point system;
% A:     n x n matrix or linear operator, (1,1) block in the saddle-
%        point matrix;
% C:     m x m matrix (m <= n), -C is the (2,2) block of the saddle-point
%        matrix;
% M:     operator, the action of the constraint preconditioner on a
%        vector;
% opts:  [optional] struct variable with the following (possible)
%        fields:
%        atol  - absolute tolerance for CG stopping criterion
%                [default 1e-6],
%        rtol  - relative tolerance for CG stopping criterion
%                [default 1e-6],
%        itmax - maximum number of CG iterations [default n],
%        print - display info about CG iterations [default true].
%
% OUTPUT ARGUMENTS:
% x:     n-vector, first n entries of the solution;
% y:     m-vector, last m entries of the solution;
% flag:  struct variable with the following fields:
%        niters - number of CG iterations performed,
%        solved - true if residNorm <= stopTol, false otherwise (itmax
%                 attained);
% stats: struct variable with the following fields:
%        residHistory - history of 2-norm of residuals.
%
%======================================================================

    % Set problem sizes and optional arguments.
    n = size(A,1);
    m = size(C,1);
    atol = 1.0e-6;
    rtol = 1.0e-6;
    itmax = n;
    display_info = true;

    if nargin > 4
    if isfield(opts, 'atol')
      atol = opts.atol;
    end
    if isfield(opts, 'rtol')
      rtol = opts.rtol;
    end
    if isfield(opts, 'itmax')
      itmax = opts.itmax;
    end
    if isfield(opts, 'print')
      display_info = opts.print;
    end
    end

    % Initialize some vectors.
    zeron = zeros(n,1);
    zerom = zeros(m,1);
    x = zeron;
    a = zerom;
    w = zerom;
    g = -b;

    % Set initial residual norm and stop tolerance.
    ru = M * [g ; w]; r = ru(1:n); u = ru(n+1:n+m);
    p = -r;
    q = -u;
    
    residNorm2 = g' * r;
    residNorm = sqrt(residNorm2);
    stopTol = atol + rtol * residNorm;
    residHistory = [residNorm];

    itn = 0;    % iteration index

    % Print initial iteration and residual norm (if required).
    if display_info
        fprintf('\n**** Constraint-preconditioned version of CG ****\n\n');
        fprintf('%5s  %9s  %9s  %9s  %9s\n', ...
                'iter', 'resid', 'pr-curv', 'du-curv', 'steplen');
        fprintf('%5d  %9.2e  ', itn, residNorm);
    end

    % Main loop.
    while residNorm > stopTol && itn < itmax

        itn = itn + 1;

        Ap = A * p; pAp = p' * Ap;
        Cq = C * q; qCq = q' * Cq;

        alpha = residNorm2 / (pAp + qCq);

        % Print curvatures and steplength (if required).
        if display_info
            fprintf('%9.2e  %9.2e  %9.2e\n', pAp, qCq, alpha);
        end

        x = x + alpha *  p;
        a = a + alpha *  q;
        g = g + alpha * Ap;
        w = w + alpha * Cq;
        
        ru = M * [g ; w]; r = ru(1:n); u = ru(n+1:n+m);
        t = a + u;
        residNorm2_new = g' * r + t' * w;
        beta = residNorm2_new / residNorm2;

        p = -r + beta * p;
        q = -t + beta * q;
 
        residNorm2 = residNorm2_new;
        residNorm = sqrt(residNorm2);
        residHistory = [residHistory; residNorm];
        
        % Print current iteration and residual norm (if required).
        if display_info
            fprintf('%5d  %9.2e  ', itn, residNorm);
        end
        
    end
    
    if display_info
        fprintf('\n'); 
    end

    stats.residHistory = residHistory;

    % Wrap up.
    flags.niters = itn;
    flags.solved = residNorm <= stopTol;
    y = a;

end
