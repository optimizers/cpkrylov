function [x, y, stats, flag] = cpcglanczos(b, A, C, M, opts)

%======================================================================
% [x, y, stats, flag] = cpcglanczos(b, A, C, M, opts)
%
% Constraint-preconditioned Lanczos version of CG (CGLanczos) for
% regularized saddle-point systems.
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
% where A is n x n, B is m x n, C is m x m, and A and C are symmetric.
% The system must also be second-order sufficient, i.e., it must satisfy
% the following condition:
%
% Let C = EDE' be a decomposition of C with D nonsingular;
% the block-diagonal matrix blkdiag(A, inv(D)) must be positive
% definite on the nullspace of [B E].
%
% The method uses a constraint preconditioner of the form
%
%  [ G   B' ]
%  [ B  -C  ],
%
% where G is a symmetric approximation to A, and must be chosen so that
% the preconditioner is second-order sufficient, i.e., blkdiag(G, inv(D))
% is positive definite on the nullspace of [B E].
%
% The contraint-preconditioned CG iterations are the same, in exact
% arithmetic, as the standard CG iterations on a reduced system.
% Below, the operator of that reduced system is denoted `op`. It can be
% viewed as the restriction of blkdiag(A, inv(D)) to the nullspace of
% [B E].
%
% The iterations stop when one of the following conditions is satisfied:
%
%  1. |r| <= stopTol  = atol + rtol * |r0|,  or
%  2. |r| <= bstopTol = btol * (|op| * |x| + |b|), or
%  3. k == itmax,
%
% where |r| and |r0| are the 2-norm of the current and initial residuals,
% atol and rtol are absolute and relative tolerances,
% btol is a relative tolerance used to stop based on a measure of the
% backward error, |op| and |x| are estimated along the iterations,
% k is the iteration index, and itmax is the maximum number of iterations.
%
% By default, btol = 0, i.e., only stopping conditions 1 and 3 are used.
%
% NOTE that B is not explicitly passed to cpcglanczos as an argument, but it
% will have been used to form the constraint preconditioner stored in M
% (see reg_cpkrylov.m).
%
% Linear operators are defined using the Spot Toolbox by Ewout van
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
% A:     n x n matrix or linear operator;
% C:     m x m matrix or linear operator;
% M:     matrix or operator materializing the action of the constraint
%        preconditioner on a vector, i.e., M*z returns the solution of
%          [ G   B' ] u = z;
%          [ B  -C  ]
% opts:  [optional] struct variable with the following fields:
%        atol  - absolute tolerance for stopping criterion 1
%                [default 1e-6],
%        rtol  - relative tolerance for stopping criterion 1
%                [default 1e-6],
%        btol  - relative tolerance used in stopping criterion 2
%                (backward error) [default 0],
%        itmax - maximum number of CP-CGLanczos iterations [default n],
%        print - display info about CP-CGLanczos iterations [default true].
%
% OUTPUT ARGUMENTS
% x:     n-vector, first n entries of the solution;
% y:     m-vector, last m entries of the solution;
% stats: struct variable with the following fields:
%        niters - number of CP-CGLanczos iterations performed,
%        residHistory - (niters+1)-vector, history of 2-norm of
%                 residuals,
%        status - string, reason why the algorithm stopped, i.e.
%                 'residual small compared to initial residual', or
%                 'backward error small', or
%                 'maximum number of iterations';
% flag:  struct variable with the following fields:
%        solved - true if |r| <= stopTol, or |r| <= bstopTol,
%                 false otherwise (itmax attained).
%======================================================================

    % Set problem sizes and optional arguments.
    n = size(A, 1);
    m = size(C, 1);
    atol = 1.0e-6;
    rtol = 1.0e-6;
    btol = 0.0;
    itmax = n;
    display_info = true;

    if nargin > 4
        if isfield(opts, 'atol')
            atol = opts.atol;
        end
        if isfield(opts, 'rtol')
            rtol = opts.rtol;
        end
        if isfield(opts, 'btol')
            btol = opts.btol;
        end
        if isfield(opts, 'itmax')
            itmax = opts.itmax;
        end
        if isfield(opts, 'print')
            display_info = opts.print;
        end
    end

    % Set up zero vectors.
    zeron = zeros(n, 1);
    zerom = zeros(m, 1);

    % Initialize some vectors, including (fake) Lanczos vectors.
    % v0 and q0.
    x = zeron;
    y = zerom;                % yk = y0 - qk, y0 = 0 ==> yk = -qk
    u = b;                    % u0 = b - A*x0 = b
    t = zerom;                % t0 = C * q0 = 0
    vk = zeron;
    qk = zerom;
    oldbeta = 0;              % used to estimate the norm of the reduced operator
    opNorm2 = 0;              % estimate of the norm of the reduced operator

    if display_info
        fprintf('\n**** Constraint-preconditioned version of CP-CGLanczos ****\n\n');
    end

    % Set Lanczos vectors v1 and q1, and initial residual norm.
    vprec = M * [u; t];       % M * [u ; -t0], t0 = 0
    vkp1  = vprec(1:n);
    qkp1  = - vprec(n+1:n+m); % q1   = q0 - vprec(n+1:n+m) = - vprec(n+1:n+m)
    beta  = dot(u, vkp1);     % beta = dot(u0, v1) + dot(t0, q1), t0 = 0
    eps100 = 100*eps;
    if beta < -eps100
        betastr = num2str(beta);
        exc = MException('CPCGLanczos:IndefiniteError', ...
                         sprintf(['Iter 0, beta (before sqrt) = ' betastr ' : preconditioner not second-order sufficient']));
        throw(exc);
    else
        bbeta = abs(beta);
        if bbeta >= eps100
            % Normalize Lanczos vectors v1 and q1.
            beta = sqrt(bbeta);
            vkp1 = vkp1 / beta;
            qkp1 = qkp1 / beta;
        end
    end
    wv = vkp1;
    wq = qkp1;
    beta1 = beta;
    residNorm = beta1;
    residHistory = [residNorm];

    % For backward error stopping criterion.
    beta1 = beta;

    % Misc. initializations.
    k    = 0;                 % iteration index
    dg   = 0;                 % d0
    low  = 1;                 % l1
    eta  = beta;              % eta1

    % Quantities related to norm(x).
    rhobar = 1;
    xxNorm2 = 0;
    xNorm = 0;
    tau = 0;
    delta = 0;

    % Stopping criterion based on residual.
    stopTol = atol + rtol * residNorm;

    % Backward error stopping criterion.
    bstopTol = btol * beta1;

    % Print initial iteration and residual norm (if required).
    if display_info
        header_fmt = '%5s  %9s';
        info_fmt = '%5d  %9.2e';
        fprintf('stopTol = %e, bstopTol = %e\n',stopTol, bstopTol);
        fprintf(header_fmt, 'iter', '|resid|');
        if btol > 0
            % Make room to print backward error, operator norm and iterate norm.
            header_fmt2 = '  %9s  %9s  %9s';
            fprintf(header_fmt2, 'bkerr', '|op|', '|x|');
        end
        fprintf('\n');
        fprintf(info_fmt, k, residNorm);
        if btol > 0
            info_fmt2 = '  %9.2e  %9.2e';
            fprintf(info_fmt2, 0, 0, xNorm);
        end
        fprintf('\n');
    end

    % Main loop.
    while residNorm > stopTol && residNorm > bstopTol && k < itmax

        k = k + 1;

        % Shift position of Lanczos vectors.
        vkm1 = vk;
        qkm1 = qk;
        vk = vkp1;
        qk = qkp1;

        % Update x and y.
        u = A * vk;
        t = C * qk;
        alpha = dot(u, vk) + dot(t, qk);

        dg = alpha - low * low * dg;      % dk
        zeta = eta / dg;                  % zetak
        x = x + zeta * wv;
        y = y - zeta * wq;                % qk = qk-1 + zetak*wqk, yk = y0-qk = -qk

        % Compute next Lanczos vectors and update residual norm.
        vprec = M * [u; -t];
        vkp1 = vprec(1:n) - alpha * vk - beta * vkm1;
        qkp1 = qk - vprec(n+1:n+m);
        qkp1 = qkp1 - alpha * qk - beta * qkm1;
        bb1 = dot(u, vkp1); bb2 = dot(t, qkp1);
        beta = dot(u, vkp1) + dot(t, qkp1);
        if beta < -eps100
            betastr = num2str(beta);
            itstr = num2str(k);
            errmsg = ['Iter ' itstr ', beta (before sqrt) = ' betastr ' : preconditioner not second-order sufficient'];
            exc = MException('CPCGLanczos:IndefiniteError', ...
                         sprintf(errmsg));
            throw(exc);
        else
            bbeta = abs(beta);
            if bbeta >= eps100
                beta = sqrt(bbeta);
                vkp1 = vkp1 / beta;
                qkp1 = qkp1 / beta;
            end
        end

        % Compute data for next updates of x and y.
        low = beta / dg;                     % lk+1
        eta = -low * eta;                    % etak+1
        wv = vkp1 - low * wv;                % wvk+1
        wq = qkp1 - low * wq;                % wqk+1

        % Compute norm(x) if using backward error stopping condition.
        if btol > 0
            rho = sqrt(rhobar * rhobar + low * low);
            cs = rhobar / rho;
            sn = low / rho;

            num = zeta - delta * tau;
            taubar = num / rhobar;
            tau = num / rho;

            xNorm = sqrt(xxNorm2 + taubar * taubar);
            xxNorm2 = xxNorm2 + tau * tau;       % for the next iteration
            delta = sn;
            rhobar = -cs;

            % Estimate norm of reduced operator.
            opNorm2 = opNorm2 + alpha * alpha + beta * beta + oldbeta * oldbeta;
            opNorm  = sqrt(opNorm2);

            bkerr = opNorm * xNorm + beta1;
            bstopTol = btol * bkerr;
        end

        residNorm = beta * abs(zeta);
        residHistory = [residHistory; residNorm];
        oldbeta = beta;

        % Print current iteration and residual norm (if required).
        if display_info
            fprintf(info_fmt, k, residNorm);
            if btol > 0
                fprintf(info_fmt2, residNorm / bkerr, opNorm, xNorm);
            end
            fprintf('\n');
        end
    end

    if display_info
        fprintf('\n');
    end

    % Wrap up.
    stats.niters = k;
    stats.residHistory = residHistory;
    flag.solved = false;
    stats.status = 'maximum number of iterations attained';
    if residNorm <= stopTol
            flag.solved = true;
            stats.status = 'residual small compared to initial residual';
    end
    if btol > 0
        if residNorm <= bstopTol
            flag.solved = true;
            stats.status = 'backward error small';
        end
    end
end
