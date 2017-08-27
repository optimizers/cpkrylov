function [x, y, flags, stats] = cpcg(A, B, C, b, opts)
  %
  % [x, resid, itn] = cpcg(H, C, P, g, maxit, atol, rtol)
  %
  % Constraint-preconditioned CG.
  % Solve the regularized saddle-point system
  %
  %  [ A   B' ] [x] = [b]
  %  [ B  -C  ] [y]   [0]
  %
  % satisfying the following conditions:
  %
  % Suppose C = EDE' is a decomposition of C with D nonsingular.
  % The block-diagonal matrix blkdiag(H, inv(D)) must be positive
  % definite on the nullspace of [A E].
  %
  % The method uses a constraint preconditioner of the form
  %
  %  [ G   A' ]
  %  [ A  -C  ],
  %
  % where G is an approximation of H, and must be chosen so that
  % blkdiag(G, inv(D)) is positive definite on the nullspace of [A E].
  %
  % The input argument P must be an operator such that P*z returns the
  % solution of
  %
  %  [ G   A' ] [r] = [z1]
  %  [ A  -C  ] [u]   [z2].
  %
  % See:
  % H.S. Dollar, N.I.M. Gould, W.H.A. Schilders, and A.J. Wathen,
  % Implicit-Factorization Preconditioning and Iterative Solvers for
  % Regularized Saddle-Point Systems, SIAM Journal on Matrix Analysis
  % and Applications, 28(1), pp. 170-189, 2006.
  %

  % Set problem sizes and optional arguments.
  n = size(A,1);
  m = size(C,1);
  atol = 1.0e-8;
  rtol = 1.0e-6;
  itmax = 2*n;         % To be modified?
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

  x = zeros(n,1);
  a = zeros(m,1);
  w = zeros(m,1);
  g = -b(1:n,1);

  ru = M * [g ; w]; r = ru(1:n); u = ru(n+1:n+m);
  p = -r;
  q = -u;
  residNorm = g' * r;
  stopTol = atol + rtol * residNorm;
  residHistory = [residNorm];

  itn = 0;

  if display_info
    fprintf('\n');
    fprintf('%4s  %8s  %8s  %8s  %8s\n', ...
            'iter', 'resid', 'pr-curv', 'du-curv', 'step');
    fprintf('%4d  %8.1e  ', itn, residNorm);
  end

  while residNorm > stopTol && itn < itmax

    itn = itn + 1;

    Ap = A * p; pAp = p' * Ap;
    Cq = C * q; qCq = q' * Cq;

    alpha = residNorm / (pAp + qCq);

    if display_info
      fprintf('%8.1e  %8.1e  %8.1e\n', pAp, qCq, alpha);
    end

    x = x + alpha *  p;
    a = a + alpha *  q;
    g = g + alpha * Ap;
    w = w + alpha * Cq;

    ru = M * [g ; w]; r = ru(1:n); u = ru(n+1:n+m);

    t = a + u;
    residNorm_new = g' * r + t' * w;
    beta = residNorm_new / residNorm;

    p = -r + beta * p;
    q = -t + beta * q;

    residNorm = residNorm_new;
    residHistory = [residHistory; residNorm];
    if display_info, fprintf('%4d  %8.1e  ', itn, residNorm); end
  end
  if display_info, fprintf('\n'); end

  stats.residHistory = residHistory;

  % Wrap up.
  flags.niters = itn;
  flags.solved = residNorm <= stopTol;
  y = -q;
end
