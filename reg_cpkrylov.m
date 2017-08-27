function [x, flags, stats] = reg_cpkrylov(method, b, A, B, C, G, opts)
  %
  % [x, flags, stats] = reg_cpkrylov(method, b, A, B, C, G, opts)
  %
  % Run a constraint-preconditioned variant of `method` for the solution of
  % block linear systems of the form
  %
  %    [ A  B' ] [ x1 ] = [ b1 ]
  %    [ B  -C ] [ x2 ]   [ b2 ]
  %
  % where A is square but need not be symmetric and C is square and
  % symmetric. The method uses a preconditioner of the form
  %
  %    [ G  B' ]
  %    [ B  -C ]
  %
  % where G is symmetric. The argument A may be a linear operator but C, B
  % and G must be explicit matrices. Let C = EFE' where F is square of
  % size rank(C). The matrix
  %
  %    [ G    0   ]
  %    [ 0 inv(F) ]
  %
  % should also be positive definite on the nullspace of J = [B E], i.e.,
  % z'Gz > 0 for all nonzero z in the nullspace of J.
  %
  % The argument `method` should be a function handle.
  %
  % The only fields of `opts` that are recognized by `cpkrylov` are:
  %  * `nitref`    : the maximum number of iterative refinement steps
  %  * `itref_tol` : the tolerance that triggers iterative refinement
  %  * `residual_update` : perform Gould-Hribar-Nocedal residual update.
  %
  % See the documentation of `method` for a description of the other input
  % and output parameters.
  %
  % Example: reg_cpkrylov(@reg_gmres3, b, A, B, C, speye(size(A)));
  %
  % Reference:
  % D. Di Serafino and D. Orban,
  % Regularized Constraint-Preconditioned Krylov Methods for General Saddle-Point Systems.
  % TBA
  %
  % dominique.orban@gerad.ca, 2017
  % daniela.diserafino@unina2.it, 2017

  % Set up coefficient matrix and constraint preconditioner.
  n = size(A,1);
  m = size(B,1);
  LDL = opLDL2(G, B, -C);

  % Check optional arguments.
  if nargin > 6
    if isfield(opts, 'nitref')
      LDL.nitref = opts.nitref;
    end
    if isfield(opts, 'itref_tol')
      LDL.itref_tol = opts.itref_tol;
    end
    if isfield(opts, 'residual_update')
      LDL.residual_update = opts.residual_update;
    end
    if isfield(opts, 'force_itref')
      LDL.force_itref = opts.force_itref;
    end
  end

  opts.M = LDL;

  % Shift linear system so rhs has the form [b ; 0] and then solve it
  xy0 = LDL * [zeros(n,1); b(n+1:n+m)];
  b1 = b(1:n) - A * xy0(1:n) - B' * xy0(n+1:n+m);
  [dx, ~, flags, stats] = method(A, B, C, [b1 ; zeros(m,1)], opts);

  % Recover solution of initial system.
  x1 = xy0(1:n) + dx;
  xy = LDL * [b1 - A * dx + G * dx; zeros(m,1)];
  x2 = xy0(n+1:n+m) + xy(n+1:n+m);
  x  = [x1; x2];
