function [x, y, stats, flag] = cpgmres(b, A, C, M, opts)

%======================================================================
% [x, y, stats, flag] = cpgmres(b, A, C, M, opts)
%
% Constraint-preconditioned GMRES(l) for regularized saddle-point
% systems.
%
%======================================================================
% Last update, August 21, 2019.
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
%   (residNorm <= stopTol = atol + rtol * residNorm0)  or  (k = itmax),
%
% where residNorm and residNorm0 are the 2-norms of the current and
% initial residuals, atol and rtol are absolute and relative tolerances,
% k is the iteration index, and itmax is the maximum number of iterations.
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
%   D. di Serafino and D. Orban,
%   Constraint-Preconditioned Krylov Solvers for Regularized
%   Saddle-Point Systems,
%   Cahier du GERAD G-2019-72, GERAD, Montreal, October 2019.
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
% opts:  [optional] struct variable with the following fields:
%        atol    - absolute tolerance for CP-GMRES stopping criterion
%                  [default 1e-6],
%        rtol    - relative tolerance for CP-GMRES stopping criterion
%                  [default 1e-6],
%        restart - restart parameter l of CP-GMRES(l) [default 50],
%        reorth  - partial reorthogonalization, true/false [default false],
%                  not implemented yet
%        itmax   - maximum number of CP-GMRES iterations [default n+m],
%        print   - display info about CP-GMRES iterations [default true].
%
% OUTPUT ARGUMENTS
% x:     n-vector, first n entries of the solution;
% y:     m-vector, last m entries of the solution;
% stats: struct variable with the following fields:
%        niters - number of CP-GMRES iterations performed,
%        residHistory - history of 2-norm of residuals;
% flag:  struct variable with the following fields (for now):
%        solved - true if residNorm <= stopTol, false otherwise (itmax
%                 attained).
%
%======================================================================

    % Set problem sizes and optional arguments.
    n = size(A,1);
    m = size(C,1);
    atol = 1.0e-6;
    rtol = 1.0e-6;
    restart = 50;
    % reorth = false;
    itmax = n+m;
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
        %   if isfield(opts, 'reorth')
        %     reorth = opts.reorth;
        %   end
        if isfield(opts, 'restart')
            restart = opts.restart;
        end
        if isfield(opts, 'print')
            display_info = opts.print;
        end
    end

    % Set up workspace.
    g = zeros(restart+1, 1);
    V = zeros(n, restart+1);        % Preconditioned Krylov vectors [v1 v2 ... vk].
    Q = zeros(m, restart+1);        % Preconditioned Krylov vectors [q1 q2 ... qk].
    H = zeros(restart+1, restart);  % Upper Hessenberg form of the saddle-point matrix.
    c = zeros(restart, 1);          % Givens cosines.
    s = zeros(restart, 1);          % Givens sines.

    % Set up zero vectors
    zeron = zeros(n,1);
    zerom = zeros(m,1);

    % Set starting guess
    x = zeron;
    y = zerom;

    % Misc. initializations.
    finished = false;
    outer = 0;
    outermax = ceil(itmax/restart);
    if display_info
        fprintf('\n**** Constraint-preconditioned version of GMRES(%d) ****\n\n', restart);
        header_fmt = '%5s  %9s\n';
        %info_fmt = '%5d  %9.2e\n';
        info_fmt = '%5d  %14.7e\n';
    end
    
    % Outer loop.
    while ~finished && outer < outermax

        outer = outer + 1;

        % Set initial Krylov vector and residual norm.
        % Separating the case outer = 1 saves some computation. Is it
        % worthwhile?
        q = zerom;
        if outer == 1
            u = b;                % u_0 = b - A * x_0 = b
            t = zerom;            % t_0 =     C * q0 = 0, y0 = q0 = 0
            w = M * [u; -t];
            V(:,1) = w(1:n,1);
            Q(:,1) = - w(n+1:n+m,1);
        else
            u = b - A*x;
            t = C*y;
            w = M * [u; -t];
            V(:,1) = w(1:n,1);
            Q(:,1) = y - w(n+1:n+m,1);
        end
        residNorm = sqrt(dot(u, V(:,1)) + dot(t, Q(:,1)));
        if residNorm ~= 0
            V(:,1) = V(:,1) / residNorm;
            Q(:,1) = Q(:,1) / residNorm;
        end
        if outer == 1
            stopTol = atol + rtol * residNorm;
            residHistory = [residNorm];
        end

        k = 0;             % Inner iteration index.
        g(1) = residNorm;

        % Print initial iteration and residual norm (if required).
        if display_info
            printf('stopTol = %e\n',stopTol);
            if outer == 1
                fprintf(header_fmt, 'iter', '|resid|');
            end
            fprintf(info_fmt, k, residNorm);
        end

        % Inner loop.
        while residNorm > stopTol && k < restart

            k = k + 1;

            % Compute next Krylov vectors from modified Gram-Schmidt
            % process.
            u = A * V(:,k);
            t = C * Q(:,k);
            w = M * [u; -t];
            V(:,k+1) = w(1:n,1);
            Q(:,k+1) = Q(:,k) - w(n+1:n+m,1);
            for j = 1 : k
                H(j,k) = dot(V(:,j), u) + dot(Q(:,j), t);
                V(:,k+1) = V(:,k+1) - H(j,k) * V(:,j);
                Q(:,k+1) = Q(:,k+1) - H(j,k) * Q(:,j);
            end
            H(k+1,k) = sqrt(dot(u, V(:,k+1)) + dot(t, Q(:,k+1)));

            if H(k+1,k) ~= 0          % Lucky breakdown if = 0.
                V(:,k+1) = V(:,k+1) / H(k+1,k);
                Q(:,k+1) = Q(:,k+1) / H(k+1,k);
            end

            % Apply previous (symmetric) Givens rotations.
            for j = 1 : k-1
                Hjk = c(j) * H(j,k) + s(j) * H(j+1,k);
                H(j+1,k) = s(j) * H(j,k) - c(j) * H(j+1,k);
                H(j,k) = Hjk;
            end
            
            % Compute and apply current (symmetric) Givens rotation:
            % [ck  sk] [H(k,k)  ] = [*]
            % [sk -ck] [H(k+1,k)]   [0]
            [c(k), s(k), H(k,k)] = SymGivens(H(k,k), H(k+1,k));
            H(k+1,k) = 0;
            g(k+1) = s(k) * g(k);
            g(k)   = c(k) * g(k);
            residNorm = abs(g(k+1));
            residHistory = [residHistory; residNorm];

            % Print current iteration and residual norm (if required).
            if display_info
                fprintf(info_fmt, (outer - 1) * restart + k, residNorm);
            end
        end

        % Update x and y with solution of upper triangular system.
        z = H(1:k,1:k) \ g(1:k);
        x = x + V(:,1:k) * z;
        q = q + Q(:,1:k) * z;
        y = y - q;

        finished = residNorm <= stopTol;
        
    end

    % Wrap up.
    stats.niters = (outer - 1) * restart + k;
    stats.residHistory = residHistory;
    flag.solved = (residNorm <= stopTol);

end
