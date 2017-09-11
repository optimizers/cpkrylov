function [x, y, flags, stats] = cpdqgmres(b, A, C, M, opts)

%======================================================================
% [x, y, flags, stats] = cpdqgmres(b, A, C, M, opts)
%
% Constraint-preconditioned DQGMRES for generalized saddle-point
% systems.
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
% where A is n x n, B is m x n, and C is symmetric of size m x m, with
% m <= n. A need not be symmetric. The saddle-point matrix is assumed
% to be nonsingular.
%
% The method uses a constraint preconditioner of the form
%
%  [ G   B' ]
%  [ B  -C  ],
%
% where G is a symmetric approximation to A, and must satisfy the
% following condition:
%
% let C = EDE' be a decomposition of C with D nonsingular;
% the block diagonal matrix blkdiag(G, inv(D)) is positive definite
% on the nullspace of [B E].
%
% The iterations stop when
%
%   sqrt(max(1, k-mem+1)) * residNorm <= stopTol = atol + rtol * residNorm0
% or
%   k = itmax,
%
% where residNorm0 is the 2-norm of the initial residual, residNorm
% is an estimate of the 2-norm of the current residual, atol and rtol
% are absolute and relative tolerances, k is the iteration index, and
% itmax is the maximum number of iterations.
% For details on the stopping criterion, see
%   Saad & Wu, DQGMRES: a Direct Quasi-minimal Residual Algorithm
%   Based on Incomplete Orthogonalization, Numerical Linear Algebra with
%   Applications, 3(4), pp. 329-343, 1996.
%
% NOTE that
% - the argument A may be a matrix or a linear operator, but C and G
%   must be explicit matrices;
% - B is not explicitly passed to cpdqgmres as an argument, but it
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
%   D. di Serafino and D. Orban,
%   Regularized Constraint-Preconditioned Krylov Solvers for General
%   Saddle-Point Systems.
%   TBA
%
%======================================================================
% INPUT ARGUMENTS
% b:     n-vector, the vector b in the rhs of the saddle-point system;
% A:     n x n matrix or linear operator, (1,1) block of the saddle-
%        point matrix;
% C:     m x m matrix (m <= n), -C is the (2,2) block of the saddle-point
%        matrix;
% M:     operator, the action of the constraint preconditioner on a
%        vector;
% opts:  [optional] struct variable with the following (possible)
%        fields:
%        atol  - absolute tolerance for DQGMRES stopping criterion
%                [default 1e-6],
%        rtol  - relative tolerance for DQGMRES stopping criterion
%                [default 1e-6],
%        mem   - memory in DQGMRES [default 10];
%        itmax - maximum number of DQGMRES iterations [default n+m],
%        print - display info about DQGMRES iterations [default true].
%
% OUTPUT ARGUMENTS
% x:     n-vector, first n entries of the solution;
% y:     m-vector, last m entries of the solution;
% flag:  struct variable with the following fields:
%        niters - number of DQGMRES iterations performed,
%        solved - true if residNorm <= stopTol, false otherwise (itmax
%                 attained);
% stats: struct variable with the following fields:
%        residHistory - history of 2-norm of residuals.
%
%======================================================================

    % Set problem sizes and optional arguments.
    n = size(A,1);
    m = size(C,1);
    atol = 1.0e-8;
    rtol = 1.0e-6;
    itmax = n+m;          % which value?????
    mem = 10;
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
        if isfield(opts, 'mem')
            mem = max(1, opts.mem);
        end
        if isfield(opts, 'print')
            display_info = opts.print;
        end
    end

    % Set up workspace.
    mem = min(mem, itmax);
    g  = zeros(mem+1, 1);
    V  = zeros(n, mem+1);     % Preconditioned Krylov vectors [v1 v2 ... vk].
    Q  = zeros(m, mem+1);     % Preconditioned Krylov vectors [q1 q2 ... qk].
    PV = zeros(n, mem+1);     % Update directions for x: PV := V * inv(R).
    PQ = zeros(m, mem+1);     % Update directions for q: PQ := Q * inv(R).
    H  = zeros(mem+1, mem);   % Upper Hessenberg form of A.
    c  = zeros(mem, 1);       % Givens cosines.
    s  = zeros(mem, 1);       % Givens sines.

    % Set up zero vectors
    zeron = zeros(n,1);
    zerom = zeros(m,1);

    % Initialize some vectors
    x = zeron;
    y = zerom;
    u = b;                    % u_0 = b - A*x_0 = b
    t = zerom;                % t_0 = C * q_0 = 0

    % Set Lanczos vectors v1 and q1, and initial residual norm.
    w = M * [u; t];           % M * [u0; -t0], t0 = 0
    V(:,1) = w(1:n,1);
    Q(:,1) = - w(n+1:n+m,1);
    residNorm = sqrt(dot(u, V(:,1))); % residnorm = sqrt(dot(u0, v1) + dot(t0, q1));
    if residNorm ~= 0
        V(:,1) = V(:,1) / residNorm;
        Q(:,1) = Q(:,1) / residNorm;
    end

    % Misc. initializations.
    k = 0;                     % Iteration index.
    g(1) = residNorm;
    stopTol = atol + rtol * residNorm;
    residHistory = [residNorm];

    % Print initial iteration and residual norm (if required).
    if display_info
        fprintf('\n**** Constraint-preconditioned version of DQGMRES - mem = %d ****\n\n', mem);
        header_fmt = '%5s  %9s\n';
        info_fmt = '%5d  %9.2e\n';
        fprintf(header_fmt, 'iter', '|resid|');
        fprintf(info_fmt, k, residNorm);
    end

    % Main loop.
    
    % The following stopping criterion compensates for the lag in the
    % residual, but usually increases the number of iterations. 
    % while sqrt(max(1, k-mem+1)) * residNorm > stopTol && k < itmax
    
    while residNorm > stopTol && k < itmax  % less accurate, but acceptable

        k = k + 1;

        % Set position in circular stack where (k+1)-st Krylov vector should go.
        kpos = mod(k-1, mem+1) + 1;  % Position corresponding to k in the circular stack.
        kp1pos = mod(k, mem+1) + 1;  % Position corresponding to k+1 in the circular stack.

        % Compute next Krylov vectors from the modified Gram-Schmidt process.
        % Only orthogonalize against the most recent min(k,mem) vectors.
        u = A * V(:,kpos);
        t = C * Q(:,kpos);
        w = M * [u; -t];
        V(:,kp1pos) = w(1:n,1);
        Q(:,kp1pos) = Q(:,kpos) - w(n+1:n+m,1);
        for j = max(1, k-mem+1) : k
            jpos = mod(j-1, mem+1) + 1;
            H(jpos,kpos) = dot(V(:,jpos), u) + dot(Q(:,jpos), t);
            V(:,kp1pos) = V(:,kp1pos) - H(jpos,kpos) * V(:,jpos);
            Q(:,kp1pos) = Q(:,kp1pos) - H(jpos,kpos) * Q(:,jpos);
        end
        H(kp1pos,kpos) = sqrt(dot(u, V(:,kp1pos)) + dot(t, Q(:,kp1pos)));

        if H(kp1pos,kpos) ~= 0     % Lucky breakdown if = 0.
            V(:,kp1pos) = V(:,kp1pos) / H(kp1pos,kpos);
            Q(:,kp1pos) = Q(:,kp1pos) / H(kp1pos,kpos);
        end

        % Apply previous (symmetric) Givens rotations.
        for j = max(1,k-mem+1) : k-1
            jpos = mod(j-1, mem+1) + 1;
            jp1pos = mod(j, mem+1) + 1;
            Hjk = c(jpos) * H(jpos,kpos) + s(jpos) * H(jp1pos,kpos);
            H(jp1pos,kpos) = s(jpos) * H(jpos,kpos) - c(jpos) * H(jp1pos,kpos);
            H(jpos,kpos) = Hjk;
        end

        % Compute and apply current (symmetric) Givens rotation:
        % [ck  sk] [H(k,k)  ] = [*]
        % [sk -ck] [H(k+1,k)]   [0].
        [c(kpos), s(kpos), H(kpos,kpos)] = SymGivens(H(kpos,kpos), H(kp1pos,kpos));
        H(kp1pos,kpos) = 0;
        g(kp1pos) = s(kpos) * g(kpos);
        g(kpos)   = c(kpos) * g(kpos);

        % Update directions PV and PQ, and solution [x; y].
        PV(:,kpos) = V(:,kpos);
        PQ(:,kpos) = Q(:,kpos);
        for j = max(1,k-mem) : k-1
            jpos = mod(j-1, mem+1) + 1;
            PV(:,kpos) = PV(:,kpos) - H(jpos,kpos) * PV(:,jpos);
            PQ(:,kpos) = PQ(:,kpos) - H(jpos,kpos) * PQ(:,jpos);
        end
        PV(:,kpos) = PV(:,kpos) / H(kpos,kpos);
        PQ(:,kpos) = PQ(:,kpos) / H(kpos,kpos);
        x = x + g(kpos) * PV(:,kpos);
        y = y - g(kpos) * PQ(:,kpos);

        % Update residual norm estimate
        residNorm = abs(g(kp1pos));
        residHistory = [residHistory; residNorm];

        % Print current iteration and residual norm (if required).
        if display_info
            fprintf(info_fmt, k, residNorm);
        end
    end

    stats.residHistory = residHistory;

    % Wrap up.
    flags.niters = k;
    flags.solved = residNorm <= stopTol;

end
