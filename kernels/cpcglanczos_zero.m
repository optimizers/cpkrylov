function [x, y, flags, stats] = cpcglanczos2(b, A, B, C, M, opts)
%======================================================================
    % [x, y, flags, stats] = cpcglanczos(b, A, C, M, opts)
    %
    % Constraint-preconditioned Lanczos version of CG (CG-Lanczos)
    %
    % Solve the regularized saddle-point system
    %
    %  [ A   B' ] [x] = [b]
    %  [ B  -C  ] [y]   [0]
    %
    % satisfying the following conditions:
    %
    % suppose C = EDE' is a decomposition of C with D nonsingular;
    % the block-diagonal matrix blkdiag(A, inv(D)) must be positive
    % definite on the nullspace of [B E].
    %
    % The method uses a constraint preconditioner of the form
    %
    %  [ G   B' ]
    %  [ B  -C  ],
    %
    % where G is an approximation of H, and must be chosen so that
    % blkdiag(G, inv(D)) is positive definite on the nullspace of [B E].
    %
    % Note that B is not explicitly passed to cpcglanczos as an argument,
    % but it has been used to form the constraint preconditioner stored
    % in M (see reg_cpkrylov.m). M must be an operator such that M*z
    % returns the solution of
    %
    %  [ G   B' ] [r] = [z1]
    %  [ B  -C  ] [u]   [z2].
    %
    % The iterations stop either when
    %
    %          residNorm <= stopTol = atol + rtol * residNorm0,
    % 
    % where residNorm and residNorm0 are the 2-norms of the current and
    % initial residuals, and atol and rtol are absolute and relative
    % tolerances, or when a maximum number of iterations is achieved.
    %
    % Reference:
    %   D. di Serafino and D. Orban,
    %   Regularized Constraint-Preconditioned Krylov Solvers for General
    %   Saddle-Point Systems.
    %   TBA
    %
    %======================================================================
    % daniela.diserafino@unicampania.it, 2017
    % dominique.orban@gerad.ca, 2017.
    %
    %======================================================================
    % INPUT ARGUMENTS
    % b:     n-vector, the vector b in the rhs of the saddle-point system;
    % A:     nxn matrix or linear operator, the (1,1) block in the saddle 
    %        point matrix;
    % C:     mxn matrix (m<=n), -C is the (2,2) block in the saddle-point
    %        matrix;
    % M:     operator, the action of the constraint preconditioner on a
    %        vector;
    % opts:  [optional] struct variable with the following (possible)
    %        entries:
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
    % flag:  struct variable with the following entries:
    %        niters - number of CG-Lanczos iterations performed,
    %        solved - true if residNorm <= stopTol, false otherwise;
    % stats: struct variable with the following entries:
    %        residHistory - history of 2-norm of residuals.
    %
    %======================================================================


% - differs from cpcglanczos in the computation of u, t, qkp1 (this
%   computation is theoretically equivalent).
%
% =========================================================================
%
% daniela.diserafino@unicampania.it, dominique.orban@gerad.ca, 2017.


% Set problem sizes and optional arguments.
n = size(A,1);
m = size(C,1);
atol = 1.0e-8;
rtol = 1.0e-6;
itmax = 2*n;                  % To be modified?
display_info = true;
M = opEye(n+m);
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
    if isfield(opts, 'M')
        M = opts.M;
    end
    if isfield(opts, 'print')
        display_info = opts.print;
    end
end

% Set up zero vectors.
zeron = zeros(n,1);
zerom = zeros(m,1);

% Set initial guess, intial residual and (fake) Lanczos vector v0.
% if isfield(opts, 'x')
%     x = opts.x;
% else
    x = zeron;
% end
% if isfield(opts, 'y')
%     y = opts.y;
% else
    y = zerom;
% end
u = b(1:n,1);             % Only the first n entries of b are required.
%q = zerom;                % q_k = - (y_k - y_0)
t = zerom;                % t_k = C * q_k = 0 for k = 0
% if ( isfield(opts, 'x') || isfield(opts, 'y') )
%     u = u - A*x + B'*q;
%     t =   - B*x + C *y;
% end
vk = zeron;
qk = zerom;

% Set Lanczos vectors v1 and q1, and initial residual norm
vprec = M * [u; t];                % t = -t = 0 at iter 0
vkp1  = vprec(1:n);
qkp1  = vprec(n+1:n+m);            % OLD: qkp1  = qk - vprec(n+1:n+m);       % qk=0, can be deleted from this formula
beta  = dot(u, vkp1);              % beta  = dot(u, vkp1) + dot(t, qkp1), but t = 0
if beta < 0
    error('Preconditioner is not positive definite.');
end
if beta ~= 0
    % Normalize Lanczos vectors v1 and q1
    beta = sqrt(beta);
    vkp1 = vkp1/beta;
    qkp1 = qkp1/beta;
end
wv = vkp1;
wq = qkp1;
residNorm = beta;
residHistory = [residNorm];

% Misc. initializations
k    = 0;
dg   = 0;      % d_0
zeta = 1;      % zeta_0
low  = 1;      % l_1
eta  = beta;   % eta_1

% Set tolerance
stopTol = atol + rtol * residNorm;

% Print initial iteration and residual norm (if required).
if display_info
    header_fmt = '%6s  %7s\n';
    info_fmt = '%5d%1s  %7.1e\n';
    fprintf(header_fmt, 'iter', '|resid|');
    fprintf(info_fmt, k, '', residNorm);
end

% Main loop.
while residNorm > stopTol && k < itmax
    
    k = k + 1;

    % Shift position of Lanczos vectors
    vkm1 = vk;
    qkm1 = qk;
    vk = vkp1;
    qk = qkp1;

    % Update x and y
    u = A * vk + B' * qk;     % OLD: u = A * vk;
                              % OLD: t = C * qk;
    alpha = dot(u,vk);        % OLD: alpha = dot(u, vk) + dot(t, qk);     % STOP if alpha <= 0?
    dg = alpha - low*low * dg;           % d_k
    zeta = eta/dg;                       % zeta_k
    x = x + zeta * wv;
%     q = q + zeta * wq;
    y = y - zeta * wq;

    % Compute next Lanczos vectors and upddate residual norm
    vprec = M * [u; t];                           % OLD: vprec = M * [u; -t];
    vkp1 = vprec(1:n) - alpha*vk - beta*vkm1;
                                                  % OLD: qkp1 = qk - vprec(n+1:n+m);
    qkp1 = vprec(n+1:n+m) - alpha*qk - beta*qkm1; % OLD: qkp1 = qkp1 - alpha*qk - beta*qkm1;
    beta = dot(u, vkp1);                          % OLD: beta = dot(u, vkp1) + dot(t, qkp1);
    if beta < 0
        error('Preconditioner is not positive definite.');
    end
    if beta ~= 0
        beta = sqrt(beta);
        vkp1 = vkp1/beta;
        qkp1 = qkp1/beta;
    end

    low = beta / dg;                      % l_k+1
    eta = -low * eta;                     % eta_k+1
    wv = vkp1 - low*wv;                   % wv_k+1
    wq = qkp1 - low*wq;                   % wq_k+1

    residNorm = beta * abs(zeta);
    residHistory = [residHistory; residNorm];
    
    
    % Print current iteration and residual norm (if required).
    if display_info
        fprintf(info_fmt, k, residNorm);
    end
    
end

%y = -q;

if display_info
    fprintf('\n');
end

stats.residHistory = residHistory;

% Wrap up
flags.niters = k;
flags.solved = (residNorm <= stopTol);

end
