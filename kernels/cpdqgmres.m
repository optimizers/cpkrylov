function [x, y, flags, stats] = reg_dqgmres3(A, B, C, b, opts)
%
% [x, y, flags, stats] = reg_dqgmres3(A, B, C, b, opts)
%
% =========================================================================
% Daniela - DOC TO BE MODIFIED, JUST REMBER THAT
% the system matrix and the preconditioner are regularized saddle-point
% ones.
%
%    System                                Preconditioner
%
%    [ A  B' ] [ x ] = [ b1 ]             [ G  B' ]
%    [ B  -C ] [ y ]   [ b2 ]             [ B  -C ]
%
% NOTE: do not keep both q and y=-q, but use only y.
%
% =========================================================================
%
% Apply preconditioned truncated GMRES to the n-by-n system Ax=b. The preconditioner
% is assumed to be symmetric and positive definite, i.e., this method is
% equivalent to applying the standard DQGMRES to the centrally-preconditioned
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
%   opts.mem   : memory parameter (default: 5)
%   opts.x     : initial guess as a *column* vector
%   opts.print : display info at each iteration (default: True)
%
% Returns:
%   flags.solved: true if successful, false if itmax was attained
%   flags.niters: number of iterations performed.
%
% daniela.diserafino@unicampania.it, dominique.orban@gerad.ca, 2017.

% Set problem sizes and optional arguments.
n = size(A,1);
m = size(C,1);
atol = 1.0e-8;
rtol = 1.0e-6;
itmax = 2*n;         % To be modified?
mem = 5;
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
  if isfield(opts, 'mem')
    mem = max(1, opts.mem);
  end
  if isfield(opts, 'print')
    display_info = opts.print;
  end
end

% Set up workspace.
mem = min(mem, itmax);
g  = zeros(mem+1, 1);
V  = zeros(n, mem+1);     % Preconditioned Krylov vectors [v1 v2 ... vk].
Q  = zeros(m, mem+1);     % Preconditioned Krylov vectors [q1 q2 ... qk].
PV = zeros(n, mem+1);     % Update directions for x: PV := V * inv(R).
PQ = zeros(m, mem+1);     % Update directions for q: PQ := Q * inv(R).
H  = zeros(mem+1, mem);   % Upper Hessenberg form of A.
c  = zeros(mem, 1);       % Givens cosines.
s  = zeros(mem, 1);       % Givens sines.

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
   t = zerom;
end

% Initial Krylov vector.
u = b(1:n,1);               % Only the first n entries of b are required
if ( isfield(opts, 'x') || isfield(opts, 'y' ) )
    u = u - A*x + B'*q;     % q = - y
    t =   - B*x + C*y;      % q = - y,  t = - Bx + Cy = - Bx - Cq
end

w = M * [u; -t];
V(:,1) = w(1:n,1);
Q(:,1) = - w(n+1:n+m,1);
residNorm = sqrt(dot(u, V(:,1)) + dot(t, Q(:,1)));
if residNorm ~= 0
    V(:,1) = V(:,1) / residNorm;
    Q(:,1) = Q(:,1) / residNorm;
end

% Misc. initializations.
k = 0;       % Iteration index.
g(1) = residNorm;
stopTol = atol + rtol * residNorm;
residHistory = [residNorm];

if display_info
  header_fmt = '%5s  %7s\n';
  info_fmt = '%5d  %7.1e\n';
  fprintf(header_fmt, 'iter', '|resid|');
  fprintf(info_fmt, k, residNorm);
end

% Main loop.
% The stopping condition compensates for the lag in the residual.
while sqrt(max(1, k-mem+1)) * residNorm > stopTol && k < itmax
  k = k + 1;

  % Set position in circular stack where (k+1)-st Krylov vector should go.
  kpos = mod(k-1, mem+1) + 1;  % Position corresponding to k in the circular stack.
  kp1pos = mod(k, mem+1) + 1;  % Position corresponding to k+1 in the circular stack.

  % Compute next Krylov vectors from the modified Gram-Schmidt process.
  % Only orthogonalize against the most recent min(k,mem) vectors.
  u = A * V(:,kpos);
  t = C * Q(:,kpos);
  w = M * [u; -t];
  V(:,kp1pos) = w(1:n,1);
  Q(:,kp1pos) = Q(:,kpos) - w(n+1:n+m,1);
  for j = max(1, k-mem+1) : k
      jpos = mod(j-1, mem+1) + 1;
      H(jpos,kpos) = dot(V(:,jpos), u) + dot(Q(:,jpos), t);
      V(:,kp1pos) = V(:,kp1pos) - H(jpos,kpos) * V(:,jpos);
      Q(:,kp1pos) = Q(:,kp1pos) - H(jpos,kpos) * Q(:,jpos);
  end
  H(kp1pos,kpos) = sqrt(dot(u, V(:,kp1pos)) + dot(t, Q(:,kp1pos)));

  if H(kp1pos,kpos) ~= 0   % Lucky breakdown if = 0.
    V(:,kp1pos) = V(:,kp1pos) / H(kp1pos,kpos);
    Q(:,kp1pos) = Q(:,kp1pos) / H(kp1pos,kpos);
  end

  % Apply previous (symmetric) Givens rotations.
  for j = max(1,k-mem+1) : k-1
    jpos = mod(j-1, mem+1) + 1;
    jp1pos = mod(j, mem+1) + 1;
    Hjk = c(jpos) * H(jpos,kpos) + s(jpos) * H(jp1pos,kpos);
    H(jp1pos,kpos) = s(jpos) * H(jpos,kpos) - c(jpos) * H(jp1pos,kpos);
    H(jpos,kpos) = Hjk;
  end

  % Compute and apply current (symmetric) Givens rotation.
  % [ck  sk] [H(k,k)  ] = [*]
  % [sk -ck] [H(k+1,k)]   [0]
  [c(kpos), s(kpos), H(kpos,kpos)] = SymGivens(H(kpos,kpos), H(kp1pos,kpos));
  H(kp1pos,kpos) = 0;
  g(kp1pos) = s(kpos) * g(kpos);
  g(kpos)   = c(kpos) * g(kpos);

  % Update directions PV and PQ, and solution [x; y].
  PV(:,kpos) = V(:,kpos);
  PQ(:,kpos) = Q(:,kpos);
  for j = max(1,k-mem) : k-1
    jpos = mod(j-1, mem+1) + 1;
    PV(:,kpos) = PV(:,kpos) - H(jpos,kpos) * PV(:,jpos);
    PQ(:,kpos) = PQ(:,kpos) - H(jpos,kpos) * PQ(:,jpos);
  end
  PV(:,kpos) = PV(:,kpos) / H(kpos,kpos);
  PQ(:,kpos) = PQ(:,kpos) / H(kpos,kpos);
  x = x + g(kpos) * PV(:,kpos);
  y = y - g(kpos) * PQ(:,kpos);

  % Update residual norm estimate
  residNorm = abs(g(kp1pos));
  residHistory = [residHistory; residNorm];

  if display_info
    % Display current info.
    fprintf(info_fmt, k, residNorm);
  end
end

stats.residHistory = residHistory;

% Wrap up.
flags.niters = k;
flags.solved = residNorm <= stopTol;
