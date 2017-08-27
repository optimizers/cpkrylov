function [x, y, flags, stats] = reg_gmres3(A, B, C, b, opts)
%
% [x, y, flags, stats] = reg_gmres3(A, B, C, b, opts)
%
% =========================================================================
% DOC TO BE MODIFIED, JUST REMBER THAT
% the system matrix and the preconditioner are regularized saddle-point
% ones.
%
%    System                                Preconditioner
%
%    [ A  B' ] [ x ] = [ b1 ]             [ G  B' ] 
%    [ B  -C ] [ y ]   [ b2 ]             [ B  -C ]
%
% NOTE: we use q instead of y; note that q=0 at each restart, while
% y=y+update (see the notations of picture 1).
%
% =========================================================================
% Apply preconditioned GMRES to the n-by-n system Ax=b. The preconditioner
% is assumed to be symmetric and positive definite, i.e., this method is
% equivalent to applying the standard GMRES to the centrally-preconditioned
% system
%          L'AL y = L'b
% where LL' = inv(M) and Ly=x.
%
% opts.M is a linear operator representing the inverse of the preconditioner.
% More precisely, the product M*v should return the solution of the system
% Ky=v where K is the preconditioner. By default, opts.M is the identity.
%
% Other optional arguments are as follows:
%   opts.atol  : absolute stopping tolerance (default: 1.0e-8)
%   opts.rtol  : relative stopping tolerance (default: 1.0e-6)
%   opts.itmax : maximum number of iterations (default: 2n)
%   opts.reorth: perform partial reorthogonalization (default: false)
%   opts.x     : initial guess as a *column* vector
%   opts.print : display info at each iteration (default: True)
%
% Returns:
%   flags.solved: true if successful, false if itmax was attained
%   flags.niter : number of iterations performed.
%
% daniela.diserafino@unicampania.it, dominique.orban@gerad.ca, 2017.

% Set problem sizes and optional arguments.
n = size(A,1);
m = size(C,1);
atol = 1.0e-8;
rtol = 1.0e-6;
% reorth = false;           % never used ???
restart = 20;
itmax = 10*ceil(n/restart); % With no restart (and exact arithmetic) it
                            % should be itmax = ceil( ((n+rank(C)-m)/restart );
                            % since rank(C) <= m, we might use 
                            % itmax = ceil(n/restart), but it could
                            % be too small when restart is used.                       
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

% Misc. initializations.
finished = false;
outer = 0;
if display_info
  header_fmt = '%6s  %7s\n';
%   info_fmt = '%5d%1s  %7.1e\n';
  info_fmt = '%5d%1s  %14.8e\n';
end

% Set up zero vectors
zeron = zeros(n,1);
zerom = zeros(m,1);

% Set initial guess 
if isfield(opts, 'x')
   x = opts.x;
else
   x = zeron;
end
if isfield(opts, 'y')
   y = opts.y;
else
   y = zerom;
   t = zerom;            % if y = -q0 = 0, then t = 0;
end

% Outer loop.
while ~finished && outer < itmax
  outer = outer + 1;

  % Initial Krylov vector.
  u = b(1:n,1);             % Only the first n entries of b are required
  % To be removed -- deal with nonzero starting guess in the the driver
  % reg_cpkrylov
  if ( (isfield(opts, 'x') || isfield(opts, 'y')) && outer == 1 ) || ( outer > 1 )
    u = u - A*x - B'*y;
    t =   - B*x + C*y;
  end

  w = M * [u; -t];
  V(:,1) = w(1:n,1);
  % Q(:,1) = q - w(n+1:n+m,1);      %WRONG: qk = - (yk - y0)
  Q(:,1) = - w(n+1:n+m,1);
  residNorm = sqrt(dot(u, V(:,1)) + dot(t, Q(:,1)));
  if residNorm ~= 0
    V(:,1) = V(:,1) / residNorm;
    Q(:,1) = Q(:,1) / residNorm;
  end
  if outer == 1
    stopTol = atol + rtol * residNorm;
    residHistory = [residNorm];
  end

  k = 0;
  g(1) = residNorm;

  if display_info
    if outer == 1
      fprintf(header_fmt, 'iter', '|resid|');
    end
    fprintf(info_fmt, k, '', residNorm);
  end

  % Inner loop.
  while residNorm > stopTol && k < restart
    k = k + 1;

    % Compute next Krylov vectors from the modified Gram-Schmidt process.
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

    status = '';
    if H(k+1,k) ~= 0   % Lucky breakdown if = 0.
      V(:,k+1) = V(:,k+1) / H(k+1,k);
      Q(:,k+1) = Q(:,k+1) / H(k+1,k);
    end

    % Apply previous (symmetric) Givens rotations.
    for j = 1 : k-1
      Hjk = c(j) * H(j,k) + s(j) * H(j+1,k);
      H(j+1,k) = s(j) * H(j,k) - c(j) * H(j+1,k);
      H(j,k) = Hjk;
    end

    % Compute and apply current (symmetric) Givens rotation.
    % [ck  sk] [H(k,k)  ] = [*]
    % [sk -ck] [H(k+1,k)]   [0]
    [c(k), s(k), H(k,k)] = SymGivens(H(k,k), H(k+1,k));
    H(k+1,k) = 0;
    g(k+1) = s(k) * g(k);
    g(k)   = c(k) * g(k);
    residNorm = abs(g(k+1));
    residHistory = [residHistory; residNorm];

    if display_info
      % Display current info.
      fprintf(info_fmt, (outer - 1) * restart + k, status, residNorm);
    end
  end

  % Update x and y with solution of upper triangular system.
  z = H(1:k,1:k) \ g(1:k);
  x = x + V(:,1:k) * z;
  y = y - Q(:,1:k) * z;     % y = y0 - q

  finished = residNorm <= stopTol;
end

stats.residHistory = residHistory;

% Wrap up.
flags.niters = (outer - 1) * restart + k;
flags.solved = (residNorm <= stopTol);
