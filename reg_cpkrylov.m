function [x, stats, flag] = reg_cpkrylov(method, b, A, B, C, G, opts)

%======================================================================
% [x, stats, flag] = reg_cpkrylov(method, b, A, B, C, G, opts)
%
% Driver for constraint-preconditioned Krylov solvers for regularized
% saddle-point systems.
%
%======================================================================
% Last update, August 21, 2019.
% Daniela di Serafino, daniela.diserafino@unicampania.it.
% Dominique Orban, dominique.orban@gerad.ca.
%
%======================================================================
% This function runs the constraint-preconditioned Krylov solver
% specified by `method` to solve the regularized saddle-point linear
% system
%
%    [ A  B' ] [ x1 ] = [ b1 ]
%    [ B  -C ] [ x2 ]   [ b2 ],
%
% where A is n x n, B is m x n, and C is symmetric of size m x m, with
% m <= n. The saddle-point matrix is assumed to be nonsingular.
% See the documentation of `method` for further details on the saddle-
% point system.
%
% The method uses a preconditioner of the form
%
%    [ G  B' ]
%    [ B  -C ]
%
% where G is a symmetric approximation to A, and must satisfy the
% following condition:
%
% let C = EDE' be a decomposition of C with D nonsingular;
% the block diagonal matrix F = blkdiag(G, inv(D)) is positive definite
% on the nullspace of [B E], i.e., z'Gz > 0 for all nonzero z in the
% nullspace of [B E].
%
% NOTE that the argument A may be a matrix or a linear operator, but C, B
% and G must be explicit matrices. The preconditioner M is a linear
% operator too, such that M*z returns the solution of
%
%  [ G   B' ] [r] = [z1]
%  [ B  -C  ] [u]   [z2].
%
% The operator M also implements iterative refinement and residual update,
% as suggested (for the case C = 0) in
%   Nicholas I. M. Gould, Mary E. Hribar, and Jorge Nocedal,
%   On the solution of equality constrained quadratic programming
%   problems quadratic programming problems arising in optimization,
%   SIAM Journal on Scientific Computing, 23(4), pp. 1376?1395, 2001.
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
% method: handle to a constraint-preconditioned Krylov solver
%         (@cpcg, @cpcglanczos, @cpminres, @cpdqgmres, @cpgmres);
% b:      n-vector, the vector b in the rhs of the saddle-point system;
% A:      n x n matrix or linear operator, (1,1) block of the saddle-
%         point matrix;
% B:      m x n matrix (m <= n), (2,1) block of the saddle-point matrix;
% C:      m x m matrix (m <= n), -C is the (2,2) block in the saddle-
%         point matrix;
% M:      operator, the action of the constraint preconditioner on a
%         vector;
% opts:   [optional] struct variable with the following (possible)
%         fields:
%         atol    - absolute tolerance for GMRES stopping criterion
%                   [default 1e-6],
%         rtol    - relative tolerance for GMRES stopping criterion
%                   [default 1e-6],
%         restart - restart parameter r of GMRES(r) [default 20],
%         reorth  - partial reorthogonalization, true/false [default false],
%         itmax   - maximum number of GMRES iterations [default n+m],
%         print   - display info about GMRES iterations [default true],
%         nitref  - maximum number of iterative refinement steps in the
%                   application of the preconditioner [default 3]
%         itref_tol - tolerance that triggers iterative refinement
%                   [default 1.0e-8],
%         force_itref - force nitref steps of iterative refinement,
%                   true/false [default false],
%         residual_update - perform Gould-Hribar-Nocedal residual update,
%                   true/false, [default false].
% NOTES
% - The only fields of opts that are recognized by reg_cpkrylov are
%   nitref, itref_tol, force_itref, and residual_update;
%   the remaining fields of opts are recognized by the Krylov solver
%   implementation specified by `method`.
% - If opts.nitref > 0, the iterative refinement is performed if
%          rNorm >= opts.itref_tol * xNorm
%   or
%          opts.force_itref,
%   where rNorm and xNorm are the 2-norms of the residual and the
%   solution resulting fom the application of the preconditioner.
%
% OUTPUT ARGUMENTS
% x:     (n+m)-vector, computed solution;
% stats: struct variable with the following fields:
%        niters - number of iterations performed by the Krylov solver,
%        residHistory - history of 2-norm of residuals,
%        ptime - time for building the preconditioner,
%        stime - time for solving the preconditioned system;
% flag:  struct variable with the following fields (for now):
%        solved - true if the residual norm satisfies the stopping
%                 condition (see the doc in `method`), false otherwise
%                 (itmax attained).
%
%======================================================================

    % Check number of input parameters
    if (nargin < 6)
        errmsg = 'reg_cpkrylov: not enough inputs';
        error(errmsg);
    end

    % Set up coefficient matrix and constraint preconditioner
    tstartp = tic;
    n = size(A,1);
    m = size(B,1);
    M = opLDL2(G, B, -C);
    ptime = toc(tstartp);

    % Check optional arguments.
    if nargin > 6
        if isfield(opts, 'nitref')
            M.nitref = opts.nitref;
        end
        if isfield(opts, 'itref_tol')
            M.itref_tol = opts.itref_tol;
        end
        if isfield(opts, 'residual_update')
            M.residual_update = opts.residual_update;
        end
        if isfield(opts, 'force_itref')
            M.force_itref = opts.force_itref;
        end
    end

    tstarts = tic;

    % Shift linear system so rhs has the form [b ; 0]
    shift = false;
    if any(b(n+1:n+m))
        shift = true;
        xy0 = M * [zeros(n,1); b(n+1:n+m)];
        b1 = b(1:n) - A * xy0(1:n) - B' * xy0(n+1:n+m);
    else
        b1 = b(1:n);
    end

    % Solve linear system
    [dx, dy, stats, flag] = method(b1, A, C, M, opts);

    % Recover solution of initial system
    if shift
        x1 = xy0(1:n) + dx;
        x2 = xy0(n+1:n+m) + dy;
    else
        x1 = dx;
        x2 = dy;
    end
    x  = [x1; x2];
    
    stime = toc(tstarts);

    stats.ptime = ptime;
    stats.stime = stime;

end
