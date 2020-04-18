function [x, y, stats, flag] = cpminres(b, A, C, M, opts)

%======================================================================
% [x, y, stats, flag] = cpminres(b, A, C, M, opts)
%
% Constraint-preconditioned MINRES for regularized saddle-point
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
% where A is n x n, B is m x n, C is m x m, with m <= n, and A and C
% are symmetric.
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
% - B is not explicitly passed to cpminres as an argument, but it
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
%        atol  - absolute tolerance for CP-MINRES stopping criterion
%                [default 1e-6],
%        rtol  - relative tolerance for CP-MINRES stopping criterion
%                [default 1e-6],
%        itmax - maximum number of CP-MINRES iterations [default n],
%        print - display info about CP-MINRES iterations [default true].
%
% OUTPUT ARGUMENTS
% x:     n-vector, first n entries of the solution;
% y:     m-vector, last m entries of the solution;
% stats: struct variable with the following fields:
%        niters - number of CP-MINRES iterations performed,
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
    u = b;                    % u0 = b - A * x0 = b
    t = zerom;                % t0 =     C * q0 = 0   
    vk = zeron;
    qk = zerom;
    
    if display_info
        fprintf('\n**** Constraint-preconditioned version of MINRES ****\n\n');
    end

    % Set Lanczos vectors v1 and q1, and initial residual norm.
    vprec = M * [u; t];       % M * [u ; -t0], t0 = 0
    vkp1  = vprec(1:n);
    qkp1  = - vprec(n+1:n+m); % q1 = q0 - vprec(n+1:n+m) = - vprec(n+1:n+m)
    beta  = dot(u, vkp1);     % beta  = dot(u0, v1) + dot(t0, q1), t0 = 0
    if beta < 0
        betastr = num2str(beta);
        errmsg = ['Iter 0, beta = ' betastr ' : preconditioner does not behave as a spd matrix.'];: preconditioner does not behave as a spd matrix.';
        error(errmsg);
    end
    if beta ~= 0
        % Normalize Lanczos vectors v1 and q1.
        beta = sqrt(beta);
        vkp1 = vkp1/beta;
        qkp1 = qkp1/beta;
    end
    wv  = vkp1;
    wq  = qkp1;
    wv2 = zeron;
    wq2 = zerom;
    residNorm = beta;
    residHistory = [residNorm];

    % Misc. initializations.
    k        = 0;              % Iteration index.
    deltabar = 0;
    epsln    = 0;
    taubar   = beta;
    cs       = -1; 
    sn       = 0;

    % Set tolerance.
    stopTol = atol + rtol * residNorm;

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
        
        % Compute next Lanczos vectors.
        u = A * vk;
        t = C * qk;
        alpha = dot(u, vk) + dot(t, qk);
        vprec = M * [u; -t];
        vkp1 = vprec(1:n) - alpha*vk - beta*vkm1;
        qkp1 = qk - vprec(n+1:n+m);
        qkp1 = qkp1 - alpha*qk - beta*qkm1;
        beta = dot(u, vkp1) + dot(t, qkp1);
        if beta < 0
            betastr = num2str(beta);
            itstr = num2str(k);
            errmsg = ['Iter ' itstr ', beta = ' betastr ' : preconditioner does not behave as a spd matrix.'];
            error(errmsg);
        end
        if beta ~= 0
            beta = sqrt(beta);
            vkp1 = vkp1/beta;
            qkp1 = qkp1/beta;
        end

        % Apply previous rotation Qk-1 to get
        %   [deltak    epslnk+1] = [cs  sn][dbark    0   ]
        %   [gammabark dbark+1 ]   [sn -cs][alphak betak+1].
        oldeps   = epsln;
        delta    = cs*deltabar + sn*alpha; % delta1    = 0        deltak
        gammabar = sn*deltabar - cs*alpha; % gammabar1 = alpha1   gammabark
        epsln    =             sn*beta;    % epsln2    = 0        epslnk+1
        deltabar =           - cs*beta;    % deltabar2 = beta2    deltabark+1

        % Compute next rotation Qk and next tau.
        gamma  = norm([gammabar beta]); % gammak
        cs     = gammabar/gamma;        % ck
        sn     = beta/gamma;            % sk
        tau    = cs*taubar ;            % tauk
        taubar = sn*taubar ;            % taubark+1
        
        % Update  x and y.
        wv1 = wv2;
        wv2 = wv;
        wq1 = wq2;
        wq2 = wq;
        wv  = (vk - oldeps*wv1 - delta*wv2)/gamma;
        wq  = (qk - oldeps*wq1 - delta*wq2)/gamma;
        x   = x + tau*wv;
        y   = y - tau*wq;

       % Set residual norm.
       residNorm = taubar;
       residHistory = [residHistory; residNorm];
        
        % Print current iteration and residual norm (if required).
        if display_info
            fprintf(info_fmt, k, residNorm);
        end
        
    end

    if display_info
        fprintf('\n');
    end

    % Wrap up.
    stats.niters = k;
    stats.residHistory = residHistory;
    flag.solved = (residNorm <= stopTol);

end
