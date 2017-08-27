classdef opLDL < opSpot
%OPLDL   Operator representing the LDL factorization of a symmetric
%        indefinite matrix with optional iterative refinement.
%
%   opLDL(A) creates an operator for multiplication by the
%   inverse of the matrix A implicitly represented by its LDL
%   factorization. Optionally, iterative refinement is performed.
%   Note that A is an explicit matrix.
%
%   The following attributes may be changed by the user:
%    * nitref : the maximum number of iterative refinement steps (3)
%    * itref_tol : iterative refinement tolerance (1.0e-8)
%    * force_itref : force iterative refinement (false)
%
%   See also ldl.
%
%   D. Orban, 2013.
%
%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Properties
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  properties( SetAccess = private )
    A             % Input matrix
    P             % Permutation matrix. TODO: Make this a vector!
    L             % Lower triangular factor
    D             % Diagonal factor
    LDL           % LDL factorization as operator
    rNorm         % Residual norm (if iterative refinement is performed)
  end

  properties( SetAccess = public )
    nitref = 3     % Default max number of iterative refinement steps
    itref_tol = 1.0e-8       % Default iterative refinement tolerance
    force_itref = false  % Force nitref steps of iterative refinement
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Methods - Public
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % opLDL. Constructor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function op = opLDL(A)
       if nargin ~= 1
          error('Invalid number of arguments.');
       end

       % Get size of input matrix
       n = size(A,1);
       if n ~= size(A,2)
          error('Input matrix must be square.');
       end

       % Construct operator
       op = op@opSpot('LDL', n, n);
       op.cflag        = ~isreal(A);
       if ~issparse(A)
         op.A          = sparse(A);
       else
         op.A          = A;
       end
       [op.L, op.D, P] = ldl(op.A);
       op.P            = opMatrix(P);
       op.LDL = op.P * inv(op.L') * inv(op.D) * inv(op.L) * op.P';
       op.sweepflag    = true;
    end % function opLDL

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function op = set.nitref(op, val)
       op.nitref = max(0, round(val));
    end

    function op = sef.itref_tol(op, val)
       op.itref_tol = max(0, val);
    end

    function op = set.force_itref(op, val)
       if val ~= false & val ~= true
          op.force_itref = false;
       else
          op.force_itref = val;
       end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % transpose
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function opOut = transpose(op)
       opOut = op;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % conj
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function opOut = conj(op)
       opOut = op.P * inv(conj(op.L)') * inv(conj(op.D)) * inv(conj(op.L)) * op.P';
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ctranpose
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function opOut = ctranspose(op)
       opOut = op.P * inv(conj(op.L)') * inv(conj(op.D)) * inv(conj(op.L)) * op.P';
    end

  end % methods - public

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Methods - protected
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  methods( Access = protected )

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % multiply
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function y = multiply(op, x, mode)
       y = op.LDL * x;
       % Perform iterative refinement if necessary / requested
       if op.nitref > 0
          r = x - op.A * y;
          rNorm = norm(r);
          xNorm = norm(x);
          nit = 0;
          while nit < op.nitref & (rNorm >= op.itref_tol * xNorm | op.force_itref)
             dy = op.LDL * r;
             y = y + dy;
             r = x - op.A * y;
             rNorm = norm(r);
             nit = nit + 1;
          end
          op.rNorm = rNorm;
       end
    end % function multiply

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % divide
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function x = divide(op, b, mode)
       x = op.A * b;
    end % function divide

  end % methods - protected

end % classdef
