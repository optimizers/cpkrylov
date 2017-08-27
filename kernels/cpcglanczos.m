function [x, y, flags, stats] = reg_cglanczos3(A, B, C, b, opts)
%
% [x, y, flags, stats] = reg_cglanczos3(A, B, C, b, opts)
%
% =========================================================================
% DOC TO BE WRITTEN, JUST REMBER THAT
% the system matrix and the preconditioner are regularized saddle-point
% ones.
%
%    System                                Preconditioner
%
%    [ A  B' ] [ x ] = [ b1 ]             [ G  B' ] 
%    [ B  -C ] [ y ]   [ b2 ]             [ B  -C ]
%
%
% INSERT SHORT DESCRIPTION of the solution of the tridiagonal system
% arising from Lanczos process (so that the meaning of the variables
% used next can be recognized)
%
% NOTES:
% - do not keep both q and y=-q, but use only y;
% - B not explicitly used if starting guess = 0.
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
qkp1  = qk - vprec(n+1:n+m);       % qk=0, can be deleted from this formula
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

% Print initial residual norm (if required)
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
    u = A * vk;
    t = C * qk;
    alpha = dot(u, vk) + dot(t, qk);     % STOP if alpha <= 0?
    dg = alpha - low*low * dg;           % d_k
    zeta = eta/dg;                       % zeta_k
    x = x + zeta * wv;
%     q = q + zeta * wq;
    y = y - zeta * wq;
    
    % Compute next Lanczos vectors and upddate residual norm
    vprec = M * [u; -t];
    vkp1 = vprec(1:n) - alpha*vk - beta*vkm1;
    qkp1 = qk - vprec(n+1:n+m);
    qkp1 = qkp1 - alpha*qk - beta*qkm1;
    beta = dot(u, vkp1) + dot(t, qkp1);
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
end

%y = -q;

if display_info
    fprintf('\n');
end

stats.residHistory = residHistory;

% Wrap up
flags.niters = k;
flags.solved = (residNorm <= stopTol);
