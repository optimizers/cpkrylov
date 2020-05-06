function [x, y, stats, flag] = cpsymmlq(b, A, C, M, opts)

%======================================================================
% [x, y, stats, flag] = cpsymmlq(b, A, C, M, opts)
%
% Constraint-preconditioned SYMMLQ for regularized saddle-point
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
% NOTE that
% - the CG residual norm is used in the stopping criterion, but the
%   histories of CG, LQ and QR residual norms are provided in output;
% - the CG residual norm is one step ahead of the other norms.
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
% opts:  [optional] struct variable with the following fields:
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
%        cgresidHistory - history of 2-norm of CG residuals,                   
%        lqresidHistory - history of 2-norm of LQ residuals,
%        qrresidHistory - history of 2-norm of MINRES residuals,
%                         see Paige & Saunders, SINUM 1975,
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

    % Initialize some vectors and the iteration index.
    x  = zeron;
    y  = zerom;                % yk = y0 - qk, y0 = 0 ==> yk = -qk  
    u  = b;                    % u0 = b - A * x0 = b
    t  = zerom;                % t0 =     C * q0 = 0
    wv = zeron;
    wq = zerom;
    k  = 0;
    
    if display_info
        fprintf('\n**** Constraint-preconditioned version of SYMMLQ ****\n\n');
    end

    % Set Lanczos vectors v1 and q1, and initial residual norm.
    vprec = M * [u; t];       % M * [u ; -t0], t0 = 0
    vkp1  = vprec(1:n);
    qkp1  = - vprec(n+1:n+m); % q1   = q0 - vprec(n+1:n+m) = - vprec(n+1:n+m)
    beta1  = dot(u, vkp1);     % beta = dot(u0, v1) + dot(t0, q1), t0 = 0
    if beta1 < 0
        errmsg = 'Iter 0: preconditioner does not behave as a spd matrix.';
        error(errmsg);
    end
    if beta1 ~= 0
        % Normalize Lanczos vectors v1 and q1.
        beta1 = sqrt(beta1);
        vkp1  = vkp1/beta1;
        qkp1  = qkp1/beta1;
    end
    cgresidNorm = beta1;
    
    % Set tolerance.
    stopTol = atol + rtol * cgresidNorm;
    
    % Initialize histories of CG, LQ and MINRES residual norms.
    if cgresidNorm <= stopTol
        lqresidNorm    = beta1;
        qrresidNorm    = beta1;
        cgresidHistory = [cgresidNorm];
        lqresidHistory = [lqresidNorm];
        qrresidHistory = [qrresidNorm];
    else
        cgresidHistory = [];
        lqresidHistory = [];
        qrresidHistory = [];
    end
    
    % Set tolerance.
    stopTol = atol + rtol * cgresidNorm;

    % Print initial iteration index and residual norms (if required).
    if display_info
        fprintf('the printed |cgresid| is one iter ahead, unless the solver\n');
        fprintf('stops at iter = 0\n\n');
        header_fmt = '%5s   %9s   %9s   %9s\n';
        info_fmt   = '%5d  %9.2e   %9.2e    %9.2e\n';
        printf('stopTol = %e\n',stopTol);
        fprintf(header_fmt, 'iter', '|cgresid|', '|lqresid|', '|qrresid|');
        if cgresidNorm <= stopTol
            fprintf(info_fmt, k, cgresidNorm, lqresidNorm, qrresidNorm);
        end
    end
    
    done = ( cgresidNorm <= stopTol );
    
    if ~done
        
        % Shift position of Lanczos vectors.
        vk      = vkp1;
        qk      = qkp1;
        
        % Compute Lanczos vectors v2 and q2.
        u     = A * vk;
        t     = C * qk;
        alpha = dot(u, vk) + dot(t, qk);
        vprec = M * [u; -t];
        vkp1  = vprec(1:n) - alpha*vk;
        qkp1  = qk - vprec(n+1:n+m);
        qkp1  = qkp1 - alpha*qk;
        beta  = dot(u, vkp1) + dot(t, qkp1);
        if beta < 0
            errmsg = ['Iter 0, second Lanczos vector: preconditioner does not behave as a spd matrix.'];
            error(errmsg);
        end
        if beta ~= 0
            beta = sqrt(beta);
            vkp1 = vkp1/beta;
            qkp1 = qkp1/beta;
        end
        
        % Initialize some quantities.
        gammabar    = alpha;
        deltabar    = beta;
        epsdelzeta  = beta1;
        epsilonzeta = 0;
        bstep       = 0;
        snprod      = 1;
        matnorm2    = alpha * alpha + beta * beta;
        
        % Main loop.
        % Only cgresidNorm is checked in the stopping criterion.
        while cgresidNorm > stopTol && k < itmax
            
            % Update estimate of matrix norm, and LQ, MINRES and CG
            % residual norms.
            matnorm = sqrt(matnorm2);
            epsmat  = matnorm * eps;
            den     = gammabar;
            if den == 0
                den = epsmat;
            end
            lqresidNorm = norm([epsdelzeta; epsilonzeta]);
            qrresidNorm = snprod * beta1;
            cgresidNorm = qrresidNorm * beta / abs(den);
            lqresidHistory = [lqresidHistory; lqresidNorm];
            qrresidHistory = [qrresidHistory; qrresidNorm];
            cgresidHistory = [cgresidHistory; cgresidNorm];
            
            % Print iteration and residual norms (if required).
            if display_info
                fprintf(info_fmt, k, cgresidNorm, lqresidNorm, qrresidNorm);
            end
            
            % Update iteration index
            k = k + 1;
            
            % Update zeta and zetabar
            zetabar = epsdelzeta / den;
            zeta    = (snprod * zetabar + bstep) / beta1;
            
            % Shift position of Lanczos vectors and save previous beta.
            vkm1 = vk;
            qkm1 = qk;
            vk = vkp1;
            qk = qkp1;
            betaold = beta;
            
            % Compute next Lanczos vectors.
            u     = A * vk;
            t     = C * qk;
            alpha = dot(u, vk) + dot(t, qk);
            vprec = M * [u; -t];
            vkp1  = vprec(1:n) - alpha*vk - beta*vkm1;
            qkp1  = qk - vprec(n+1:n+m);
            qkp1  = qkp1 - alpha*qk - beta*qkm1;
            beta  = dot(u, vkp1) + dot(t, qkp1);
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
            
            % Update square of estimate of matrix norm.
            matnorm2 = matnorm2 + alpha * alpha + beta * beta + betaold * betaold;
            
            % Compute next plane rotation.
            gamma    = norm([gammabar; betaold]);    % gammak
            cs       = gammabar / gamma;             % csk+1
            sn       = betaold / gamma;              % snk+1
            delta    = cs * deltabar + sn * alpha;   % deltak+1
            gammabar = sn * deltabar - cs * alpha;   % gammabark+1
            epsilon  = sn * beta;                    % epsilonk+2
            deltabar =               - cs * beta;    % deltabark+2
            
            % Update x and y (LQ solution).
            zeta     = epsdelzeta / gamma;
            zcs      = zeta * cs;
            zsn      = zeta * sn;
            x        = x + zcs * wv + zsn * vk;
            y        = y - zcs * wq - zsn * qk;
            wv       = sn * wv - cs * vk;
            wq       = sn * wq - cs * qk;
            
            % Accumulate bstep and snprod and update other quantities
            % for next step
            bstep       = bstep + snprod * cs * zeta;
            snprod      = snprod * sn;
            epsdelzeta  = epsilonzeta - delta * zeta;   % -epsilonk+1 * zetak-1 + deltak+2 * zetak
            epsilonzeta = - epsilon * zeta;             % -epsilonk+2 * zetak
            
        end
        
        % Compute the last SYMMLQ and MINRES residual norms
        matnorm = sqrt(matnorm2);
        epsmat  = matnorm * eps;
        den     = gammabar;
        if den == 0
            den = epsmat;
        end
        lqresidNorm = norm([epsdelzeta; epsilonzeta]);
        qrresidNorm = snprod * beta1;
        lqresidHistory = [lqresidHistory; lqresidNorm];
        qrresidHistory = [qrresidHistory; qrresidNorm];
        
        % Align CG residual norm history with LQ and MINRES
        % residual norm histories.
        cgresidHistory = [beta1; cgresidHistory];
        
        % Move to CG solution if cgresidNorm < lqresidNorm.
        if (cgresidNorm < lqresidNorm)
            zetabar = epsdelzeta / den;
            bstep  = bstep + snprod * zetabar;
            x = x + zetabar * wv;
            y = y - zetabar * wq;
        end
        
        %  Add the step along first Lanczos vector.
        vprec = M * [b; zerom];
        vk  = vprec(1:n);
        qk  = - vprec(n+1:n+m);
        bstep  = bstep / beta1;
        x   = x + bstep * vk;
        y   = y - bstep * qk;
        
        % Print last iteration index and corresponding residual norms 
        % (if required).
        if display_info
            info_fmt   = '%5d     ---      %9.2e    %9.2e\n';
            fprintf(info_fmt, k, lqresidNorm, qrresidNorm);
        end
        
        if display_info
            fprintf('\n');
        end

    end % if ~done

    % Wrap up.
    stats.niters          = k;
    stats.lqresidHistory = lqresidHistory;
    stats.qrresidHistory = qrresidHistory;
    stats.cgresidHistory = cgresidHistory;
    flag.solved          = (cgresidNorm <= stopTol);

end
