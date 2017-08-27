classdef opLDL2 < opSpot
%OPLDL   Operator representing the LDL factorization of a symmetric
%        indefinite matrix with optional iterative refinement.
%
%   opLDL2(A,B,C) creates an operator for multiplication by the
%   inverse of the block matrix
%     [ A  B' ]
%     [ B  C  ]
%   implicitly represented by its LDL factorization.
%   Optionally, iterative refinement is performed. Note that A, B
%   and C are explicit matrices.
%
%   The following attributes may be changed by the user:
%    * nitref : the maximum number of iterative refinement steps (3)
%    * itref_tol : iterative refinement tolerance (1.0e-8)
%    * force_itref : force iterative refinement (false)
%    * residual_update : use Gould-Hribar-Nocedal residual update (false)
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
      nA            % Leading  block dimension
      nC            % Trailing block dimension
   end

   properties( SetAccess = public )
      nitref = 3    % Default max number of iterative refinement steps
      itref_tol = 1.0e-8  % Default iterative refinement tolerance
      force_itref = false  % Mostly for debugging
      Aty = false   % Previous range space component
      residual_update = false
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Methods - Public
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % opDiag. Constructor
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function op = opLDL2(A, B, C)
         if nargin ~= 3
            error('Invalid number of arguments.');
         end

         % Get size of input matrix
         nA = size(A,1);
         nC = size(C,1);
         if nA ~= size(A,2) | nC ~= size(C,2)
            error('First and last arguments must be square.');
         end
         nB = size(B,1);
         mB = size(B,2);
         if mB ~= nA | nB ~= nC
            error('Incompatible dimensions.');
         end
         n = nA + nC;

         % Construct operator
         op = op@opSpot('LDL2', n, n);
         op.cflag        = ~isreal(A);
         op.A            = [A  B' ; B  C];
         [op.L, op.D, P] = ldl(op.A);
         op.P            = opMatrix(P);
         op.LDL = op.P * inv(op.L') * inv(op.D) * inv(op.L) * op.P';
         op.sweepflag    = true;
         op.nA           = nA;
         op.nC           = nC;
      end % function opLDL2

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

      function op = set.residual_update(op, val)
         if val ~= false & val ~= true
            op.residual_update = false;
         else
            op.residual_update = val;
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

      function x = double(op)
      %double  Convert operator to a double.

         % Can't do op * eye(n), but can do op \ eye(n).
         e = zeros(op.n, 1);
         x = zeros(op.n);
         for i = 1 : op.n
           e(i) = 1;
           x(:, i) = op * e;
           e(i) = 0;
         end
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
         n = op.nA;
         m = op.nC;
         if op.residual_update & op.Aty
            y = op.LDL * [x(1:n) - op.Aty ; x(n+1:op.n)];
         else
            y = op.LDL * x;
         end
         if op.residual_update
            op.Aty = op.A(1:n, n+1:n+m) * y(n+1:n+m);
         end
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
