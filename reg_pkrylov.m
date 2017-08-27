function [x, flags, stats] = reg_pkrylov(method, b, A, B, C, G, opts)
  %
  % [x, flags, stats] = pkrylov(method, b, A, B, E, F, G, opts)
  %
  % Projected Krylov method for the solution of block linear systems
  % of the form
  %
  %    [ A   B' ] [ x1 ] = [ b1 ]
  %    [ B  -C  ] [ x2 ]   [ b2 ]
  %
  % where A is square but need not be symmetric and C is square and (FOR
  % THE MOMENT?) diagonal. The system is reformulated as
  %
  %    [ A    0     B' ] [ x1 ] = [ b1 ]
  %    [ 0  inv(F)  E' ] [ w  ] = [  0 ]
  %    [ B    E     0  ] [ x2 ]   [ b2 ]
  %
  % where EFE'= C and F is square with size rank(C). The method uses a
  % preconditioner of the form
  %
  %    [ G    0     B' ]
  %    [ 0  inv(F)  E' ]
  %    [ B    E     0  ]
  %
  % where G is symmetric. The argument A may be a linear operator, but B, C
  % and G must be explicit matrices. In addition,
  %
  %    M := [ G    0    ]
  %         [ 0  inv(F) ]
  % 
  % should be positive definite on the nullspace of [B E], i.e., z'Mz > 0
  % for all nonzero z in the nullspace of [B E].
  %
  % The argument `method` should be a function handle.
  %
  % The only fields of `opts` that recognized by `regpkrylov` are:
  %  * `nitref`    : the maximum number of iterative refinement steps
  %  * `itref_tol` : the tolerance that triggers iterative refinement
  %  * `residual_update` : perform Gould-Hribar-Nocedal residual update.
  %
  % See the documentation of `method` for a description of the other input
  % and output parameters.
  %
  % Example: pkrylov(@gmres3, b, A, B, C, G);
  %
  % References:
  % N. I. M. Gould, D. Orban and T. Rees,
  % Projected Krylov Methods for General Regularized Saddle-Point Systems.
  % ????
  %
  % dominique.orban@gerad.ca, 2017.
  % daniela.diserafino@unina2.it, 2017.
  %
  % Set up projection operator: [u1] = Proj * [w1] means
  %                             [u2]          [ z]
  %
  % [w1]    [w1]    [u1]           [ G    0     B' ] [u1]   [w1]    [u1]
  % [ z] -> [ z] -> [u2] such that [ 0  inv(F)  E' ] [u2] = [ z] -> [u2].
  %         [ 0]    [u3]           [ B    E     0  ] [u3]   [ 0]
    
  n = size(A,1);
  m = size(B,1);
  
  % B, E and C must be matrices, because GF and BE must be matrices.
  ind = find(diag(C)~=0);
  s = length(ind);
  Im = speye(m);                          % Im = opEye(m);
  E = Im(1:m,ind);                        % E = opRestriction(m,ind);
  F = C(ind,ind);                         % F = E*C*E'; C diag (for now) 
%   invF = inv(F);                          % invFop = opInverse(F);
  f = diag(F);
  invF = spdiags(1./f,0,s,s);
  AF = [A opZeros(n,s); opZeros(s,n) invF];
  GF = [G sparse(n,s); sparse(s,n) invF]; % LDL requires GF to be a matrix
  BE = [B E];
  
  % Should code LDL operator working directly on the 3x3 block form
  LDL = opLDL2(GF, BE, sparse(m,m));
  IO = [opEye(n+s,n+s) ; opZeros(m,n+s)];

  % Check optional arguments.
  if nargin > 5
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

  Proj = opRestriction(n+s+m, [1:n+s]) * LDL * IO;
  opts.M = Proj;

  % Shift linear system so rhs has the form [b1; 0 ; 0]. The second
  % subvector is already 0.
  xy0 = LDL * [zeros(n+s,1); b(n+1:n+m)];
  b1z = [b(1:n); zeros(s,1)];
  b1z = b1z - AF * xy0(1:n+s);
  [dx, flags, stats] = method(AF, b1z, opts);

  % Recover solution of initial system.
  x1z = xy0(1:n+s) + dx;
  x1 = x1z(1:n);
  xy = LDL * [b1z - AF * dx; zeros(m,1)];
  x2 = xy(n+s+1:n+s+m);
  x = [x1; x2];
