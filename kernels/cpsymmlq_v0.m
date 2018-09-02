function [x, y, stats, flag] = cpsymmlq_v0(b, A, C, M, opts)

%======================================================================
% [x, y, stats, flag] = cpsymmlq(b, A, C, M, opts)
%
% Constraint-preconditioned SYMMLQ for regularized saddle-point
% systems.
%
%======================================================================
% Last update, August 28, 2018.
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
% NOTE that
% - the CG residual norm is used in the stopping criterion, but the
%   histories of both the LQ and CG residual norms are provided in output;
% - the CG residual norm is one step ahead of the LQ residual norm
%   (see also lqresidHistory and cgresidHistory in the ouput arguments). 
%
% NOTE also that
% - the argument A may be a matrix or a linear operator, but C and G
%   must be explicit matrices;
% - B is not explicitly passed to cpsymmlq as an argument, but it
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
%        atol  - absolute tolerance for CP-SYMMLQ stopping criterion
%                [default 1e-6],
%        rtol  - relative tolerance for CP-SYMMLQ stopping criterion
%                [default 1e-6],
%        itmax - maximum number of CP-SYMMLQ iterations [default n],
%        print - display info about CP-SYMMLQ iterations [default true].
%
% OUTPUT ARGUMENTS
% x:     n-vector, first n entries of the solution;
% y:     m-vector, last m entries of the solution;
% stats: struct variable with the following fields:
%        niters         - number of CP-SYMMLQ iterations performed,
%        cgresidHistory - history of 2-norm of CG residuals
%                         (last entry corresponds to niters-th LQ iterate),
%        lqresidHistory - history of 2-norm of LQ residuals
%                         (last entry corresponds to niters-th CG iterate).                        
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

    % Initialize some vectors and iteration index
    x = zeron;
    y = zerom;                % yk = y0 - qk, y0 = 0 ==> yk = -qk  
    u = b;                    % u0 = b - A * x0 = b
    t = zerom;                % t0 =     C * q0 = 0
    k = 0;
    
    if display_info
        fprintf('\n**** Constraint-preconditioned version of SYMMLQ ****\n\n');
    end

    % Set Lanczos vectors v1 and q1, and initial residual norm.
    vprec = M * [u; t];       % M * [u ; -t0], t0 = 0
    vkp1  = vprec(1:n);
    qkp1  = - vprec(n+1:n+m); % q1   = q0 - vprec(n+1:n+m) = - vprec(n+1:n+m)
    beta  = dot(u, vkp1);     % beta = dot(u0, v1) + dot(t0, q1), t0 = 0
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
    lqresidNorm = beta;
    cgresidNorm = beta;
    lqresidHistory = [lqresidNorm];
    cgresidHistory = [cgresidNorm];

    % Set tolerance.
    stopTol = atol + rtol * cgresidNorm;

    % Print initial iteration and residual norms (if required).
    if display_info
        fprintf('Note that cgresid is one iter ahead\n\n');
        header_fmt = '%5s   %9s   %9s\n';
        info_fmt   = '%5d  %9.2e   %9.2e\n';
        fprintf(header_fmt, 'iter', '|cgresid|', '|lqresid|');
        fprintf(info_fmt, k, cgresidNorm, lqresidNorm);
    end

    if cgresidNorm > stopTol
        
        % The iteration ondex is not updated because the first rotation
        % is not applied yet and the LQ siolution is not updated
        % (possibly, only the CG solution is).

        % Shift position of Lanczos vectors and beta, and initialize wbar
        vk = vkp1;
        qk = qkp1;
        wvbar = vk;
        wqbar = qk;
        betaold = beta;
        
        % Compute Lanczos vectors v2 and q2.
        u     = A * vk;
        t     = C * qk;
        alpha = dot(u, vk) + dot(t, qk);          % STOP if alpha <= 0 ?
        vprec = M * [u; -t];
        vkp1  = vprec(1:n) - alpha*vk;
        qkp1  = qk - vprec(n+1:n+m);
        qkp1  = qkp1 - alpha*qk;
        beta  = dot(u, vkp1) + dot(t, qkp1);
        if beta < 0
            % beta
%             itstr = num2str(k);
            errmsg = ['Iter 1, second Lanczos vector: preconditioner does not behave as a spd matrix.'];
            error(errmsg);
        end
        if beta ~= 0
            beta = sqrt(beta);
            vkp1 = vkp1/beta;
            qkp1 = qkp1/beta;
        end
        
        % Compute first plane rotation.
        % gammabar and deltabar are set for uniformity with other
        % iterations, epsilonzeta is initialized for the next iteration.
        gammabar    = alpha;                     % gammabar1
        deltabar    = beta;                      % deltabar1
        gamma       = sqrt(gammabar^2 + beta^2); % gamma1
        cs          = gammabar / gamma;          % cs1
        sn          = beta / gamma;              % sn1
        zeta        = betaold / gamma;           % zeta1
        epsilonzeta = 0;                         % -epsilon2 * zeta0 (fake value)
        
        % Initialize estimate of matrix norm
        matnorm2 = alpha * alpha + beta * beta;
        matnorm  = sqrt(matnorm2);
        matnormeps = matnorm * eps;

        % Compute cgresidNorm and update CG solution if cgresidNorm
        % satisfies stopping criterion.
        cgresidNorm2 = abs(betaold * sn);
        d = cs;
        if d == 0
            d = matnormeps;
        end
        cgresidNorm  = cgresidNorm2 / abs(d);
        if cgresidNorm <= stopTol
            x = x + (zeta/cs) * wvbar;
            y = y - (zeta/cs) * wqbar; 
        end
        cgresidHistory = [cgresidHistory; cgresidNorm];

        % Print current iteration and residual norms (if required).
        if display_info
            fprintf(info_fmt, k, cgresidNorm, lqresidNorm);
        end
        
    end 
    
    % Main loop.
    % Only cgresidNorm is checked in the stopping criterion.
    while cgresidNorm > stopTol && k < itmax
        
        k = k + 1;
        
        % Update x and y (LQ solution).
        wv    = cs * wvbar + sn * vkp1;
        wq    = cs * wqbar + sn * qkp1;
        x     = x + zeta * wv;
        y     = y - zeta * wq;
        wvbar = sn * wvbar - cs * vkp1;
        wqbar = sn * wqbar - cs * qkp1;   

        % Shift position of Lanczos vectors.
        vkm1 = vk;
        qkm1 = qk;
        vk = vkp1;
        qk = qkp1;
        
        % Compute next Lanczos vectors.
        betaold = beta;
        u     = A * vk;
        t     = C * qk;
        alpha = dot(u, vk) + dot(t, qk);  % STOP if alpha <= 0 ?
        vprec = M * [u; -t];
        vkp1  = vprec(1:n) - alpha*vk - beta*vkm1;
        qkp1  = qk - vprec(n+1:n+m);
        qkp1  = qkp1 - alpha*qk - beta*qkm1;
        beta  = dot(u, vkp1) + dot(t, qkp1);
        if beta < 0
            % beta
            itstr = num2str(k);
            errmsg = ['Iter ' itstr ': preconditioner does not behave as a spd matrix.'];
            %error(errmsg);
        end
        if beta ~= 0
            beta = sqrt(beta);
            vkp1 = vkp1/beta;
            qkp1 = qkp1/beta;
        end

        % Update estimate of matrix norm.
        matnorm2 = matnorm2 + alpha * alpha + beta * beta + betaold * betaold;

        % Apply previously computed rotation and compute new one.
        % The notation follows Paige & Saunders, SINUM 1975.
        delta       =   cs * deltabar + sn * alpha;  % deltak+1
        deltazeta   = - delta * zeta;                % deltak+1 * zetak
        gammabar    =   sn * deltabar - cs * alpha;  % gammabark+1
        epsilon     =   sn * beta;                   % epsilonk+2           
        deltabar    = - cs * beta;                   % deltabark+2
        gamma       =   sqrt(gammabar^2 + beta^2);   % gammak+1
        cs          =   gammabar / gamma;            % csk+1
        sn          =   beta / gamma;                % snk+1
        epszdelz    =   epsilonzeta + deltazeta;     % -epsilonk+1 * zetak-1 + deltak+2 * zetak
        epsilonzeta = - epsilon * zeta;              % epsilonk+2 * zetak
        zeta        =   epszdelz / gamma;            % zetak+1
        
        % Compute norms of residuals corresponding to k-th LQ iterate
        % and (k+1)-st CG iterate.
        lqresidNorm  = sqrt(epszdelz^2 + epsilonzeta^2);
        cgresidNorm2 = cgresidNorm2 * abs(sn);
        d = cs;
        if d == 0
            d = matnormeps;
        end
        cgresidNorm  = cgresidNorm2 / abs(d);
        lqresidHistory = [lqresidHistory; lqresidNorm];
        cgresidHistory = [cgresidHistory; cgresidNorm];
        
        % Print current iteration and residual norms (if required).
        if display_info
            fprintf(info_fmt, k, cgresidNorm, lqresidNorm);
        end
        
    end

    if display_info
        fprintf('\n');
    end
    
    % Move to CG solution if cgresidNorm < lqresidNorm.
    % Since the convergence tests involve only cgresidNorm, it is unlikely
    % to stop at an LQ point.
    if (cgresidNorm <= lqresidNorm)
        x = x + (epszdelz/gammabar) * wvbar;
        y = y - (epszdelz/gammabar) * wqbar;
    end

    % Wrap up.
    stats.niters          = k;
    stats.cgresidHistory  = cgresidHistory;
    stats.lqresidHistory = lqresidHistory;
    flag.solved           = (cgresidNorm <= stopTol);

end
