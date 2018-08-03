function [x, y, stats, flag] = cpcglanczos2(b, A, C, M, opts)

%======================================================================
% [x, y, stats, flag] = cpcglanczos(b, A, C, M, opts)
%
% Constraint-preconditioned Lanczos version of CG (CG-Lanczos) for
% regularized saddle-point systems.
%
%======================================================================
% Last update, July 8, 2018.
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
% where G is a symmetric approximation to A, and must be chosen so that
% blkdiag(G, inv(D)) is positive definite on the nullspace of [B E].
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
% - B is not explicitly passed to cpcglanczos as an argument, but it
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
%        atol  - absolute tolerance for CG-Lanczos stopping criterion
%                [default 1e-6],
%        rtol  - relative tolerance for CG-Lanczos stopping criterion
%                [default 1e-6],
%        itmax - maximum number of CG-Lanczos iterations [default n],
%        print - display info about CG-Lanczos iterations [default true].
%
% OUTPUT ARGUMENTS
% x:     n-vector, first n entries of the solution;
% y:     m-vector, last m entries of the solution;
% stats: struct variable with the following fields:
%        niters - number of CG-Lanczos iterations performed,
%        residHistory - history of 2-norm of residuals;
% flag:  struct variable with the following fields (for now):
%        solved - true if residNorm <= stopTol, false otherwise (itmax
%                 attained).
%
%======================================================================
% NOTE
%   cpcglanczos2 differs from cpcglanczos in the update of the Lanczos
%   vector (see lines 189-201).
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

    % Set up zero vectors.
    zeron = zeros(n,1);
    zerom = zeros(m,1);

    % Initialize some vectors, including (fake) Lanczos vectors.
    % v0 and q0.
    x = zeron;
    y = zerom;                % yk = y0 - qk, y0 = 0 ==> yk = -qk
    u = b;                    % u0 = b - A*x0 = b
    t = zerom;                % t0 = C * q0 = 0
    vk = zeron;
    qk = zerom;
    
    oldbeta = 0;              % used to estimate the matrix norm
    Matnorm2 = 0;             % initialize estimate of matrix norm

    if display_info
        fprintf('\n**** Constraint-preconditioned version of CG-Lanczos ****\n\n');
    end

    % Set Lanczos vectors v1 and q1, and initial residual norm.
    vprec = M * [u; t];       % M * [u ; -t0], t0 = 0
    vkp1  = vprec(1:n);
    qkp1  = - vprec(n+1:n+m); % q1 = q0 - vprec(n+1:n+m) = - vprec(n+1:n+m)
    beta  = dot(u, vkp1);     % beta  = dot(u0, v1) + dot(t0, q1), t0 = 0
    if beta < 0
        errmsg = 'Iter 0: preconditioner does not behave as a spd matrix.';
        error(errmsg);
    end
    if beta ~= 0
        % Normalize Lanczos vectors v1 and q1.
        beta = sqrt(beta);
        vkp1 = vkp1/beta;
        qkp1 = qkp1/beta;
    end
    wv = vkp1;
    wq = qkp1;
    residNorm = beta;
    residHistory = [residNorm];
    
    % For backward error stopping criterion
    beta1 = beta;

    % Misc. initializations.
    k    = 0;                 % iteration index
    dg   = 0;                 % d0
    tau  = 0;                 % tau0
    delta = 0;                % delta0
    xnormacc = 0;             % xanormacc0
    low  = 1;                 % l1
    eta  = beta;              % eta1
    rhobar = 1;               % rhobar1
    
    % Set tolerance.
    % ADD OPTION TO CHOOSE BETWEEN STOPPING CRITERIA (RELATIVE RESIDUAL NORM
    % AND BACKWARD ERROR) AND MODIFY THIS SECTION OF CODE ACCORDINGLY
    
    % Stopping criterion based only on residual
    % stopTol = atol + rtol * residNorm;
    
    % Backward error stopping criterion
    stopTol = rtol * beta1;

    % Print initial iteration and residual norm (if required).
    if display_info
        header_fmt = '%5s  %9s\n';
        info_fmt = '%5d  %9.2e\n';
        fprintf(header_fmt, 'iter', '|resid|');
        fprintf(info_fmt, k, residNorm);
    end

    % Main loop.
    while residNorm > stopTol && k < itmax

        k = k + 1;

        % Shift position of Lanczos vectors.
        vkm1 = vk;
        qkm1 = qk;
        vk = vkp1;
        qk = qkp1;

        % Compute next Lanczos vectors, except for normalization,
        % using updating formulas of the following type
        % (improve local orthogonality in finite precision):
        % w = Mat*zk, wprec = M*w, wkp1 = wprec - beta* wkm1,
        % alpha = wkp1'*wk, wkp1 = wkp1 - alpha* zk,
        u = A * vk;
        t = C * qk;
        vprec = M * [u; -t];
        vkp1 = vprec(1:n) - beta*vkm1;
        qkp1 = qk - vprec(n+1:n+m);
        qkp1 = qkp1 - beta*qkm1;
        alpha = dot(u, vk) + dot(t, qk);  % STOP if alpha <= 0 ?
        vkp1 = vkp1 - alpha*vk;
        qkp1 = qkp1 - alpha*qk;

        % Update x and y
        dg = alpha - low*low * dg;        % dk
        zeta = eta/dg;                    % zetak
        x = x + zeta * wv;
        y = y - zeta * wq;                % qk = qk-1 + zetak*wqk, yk = y0-qk = -qk

        % Upddate residual norm and normalize Lanczos vectors
        beta = dot(u, vkp1) + dot(t, qkp1);
        if beta < 0
            itstr = num2str(k);
            errmsg = ['Iter ' itstr ': preconditioner does not behave as a spd matrix.'];
            error(errmsg);
        end
        if beta ~= 0
            beta = sqrt(beta);
            vkp1 = vkp1/beta;
            qkp1 = qkp1/beta;
        end

        % Compute data for next updates of x and y
        low = beta / dg;                     % lk+1
        eta = -low * eta;                    % etak+1
        wv = vkp1 - low*wv;                  % wvk+1
        wq = qkp1 - low*wq;                  % wqk+1
        
        % Compute data for 2-norm of [x; y]
        rho = sqrt(rhobar*rhobar + low*low); % rhok
        c = rhobar/rho;                      % ck
        s = low/rho;                         % sk
        tau = zeta - delta*tau;              % num of tauk and taubark
        taubar = tau/rhobar;                 % taubark
        tau = tau/rho;                       % tauk
        delta = s;                           % deltak+1
        rhobar = -c;                         % rhobark+1
    
        % Estimate Frobenius norm of preconditioned matrix and compute
        % 2-norm of [x; y]
        Matnorm2 = Matnorm2 + alpha * alpha + beta * beta + oldbeta * oldbeta;
        xnorm2 = xnormacc + taubar*taubar;
        xnormacc = xnormacc + tau*tau;
        xnorm = sqrt(xnorm2);
        Matnorm = sqrt(Matnorm2);
        oldbeta = beta;
        
        % Compute residual norm
        residNorm = beta * abs(zeta);
        residHistory = [residHistory; residNorm];


        % Print current iteration and residual norm (if required).
        if display_info
            info_fmt = '%5d  %9.2e  %9.2e  %16.8e  %16.8e\n';
            fprintf(info_fmt, k, residNorm, Matnorm * xnorm + beta1, xnorm, norm([x; y])); % TO BE MODIFIED
        end

        % TO BE MODIFIED: UPDATE STOPTOL ONLY IF BACKWARD ERROR STOPPING
        % CRITERION IS USED
        stopTol = rtol * (Matnorm * xnorm + beta1);
        
    end

    if display_info
        fprintf('\n');
    end

    % Wrap up.
    stats.niters = k;
    stats.residHistory = residHistory;
    flag.solved = (residNorm <= stopTol);

end

