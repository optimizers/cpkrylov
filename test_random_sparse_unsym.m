% TEST PROGRAM
% Run CPDQGMRES and CPGMRES(r) on random unsymmetric generalized
% saddle-point systems.
%
% Daniela di Serafino, Sept. 10, 2017

% clear;
% close all;

n = 200; m = 100;
density = 0.01;
Q = sprand(n,n,density,1e-1);  % Q may have zeros on diag
a = 5.0e-1;                     % the smaller a, the worse Q. Small entries on diag?
Q = Q + a*speye(n);
A = sprand(m,n,density);
% A = 1e8*A;
C = spdiags(abs(diag(rand(m))),0,m,m);
% C = 1e-2*C;
K = [Q  A' ; A  -C];
b = rand(n+m, 1);

x = K\b;

fprintf('\nSADDLE-POINT SYSTEM (n = %d, m = %d)\n', n, m);
fprintf('\n  [ A   B'' ] [x] = [b]');
fprintf('\n  [ B  -C  ] [y]   [0]\n');
fprintf('\n  matrix rank: %d\n', rank(full(K)));

% G = speye(n);
% G = spdiags(abs(diag(Q)),0,n,n);

% Define symmetric banded G with half bandwidth hbw. This does not ensure
% G is spd unless hbw = 0 (G diagonal).
hbw = 0;
qdiag = full(abs(diag(Q)));
diags = spdiags(((Q+Q'))/2, -hbw:hbw);
diags(:,hbw+1) = qdiag;
G = spdiags(diags, -hbw:hbw, n, n);
P = [G  A' ; A  -C];

% Info on the eigenvalues of the constraint-preconditioned matrix.
lambdaPK = eig(full(K),full(P));
n_lambdaPK_1 = length(find(abs(lambdaPK - 1) < 1e-8));
n_lambdaPK_inf = length(find(abs(lambdaPK) == Inf));
fprintf('\nCONSTRAINT-PRECONDITIONED MATRIX\n');
fprintf('\n  num. eigs close to 1: %d', n_lambdaPK_1);
fprintf('\n  num. eigs Inf: %d\n', n_lambdaPK_inf);
figure
plot(lambdaPK,'*');
title('Eigenvalues of preconditioned matrix');

pause

cpk3 = @cpdqgmres;
cpk_string3 = 'CPDQGMRES';
cpk4 = @cpgmres;
cpk_string4 = 'CPGMRES';

opts.print = true;     % true/false;
opts.atol = 1.0e-6;
opts.rtol = 1.0e-6;
opts.itmax = 500;
opts.mem = 10;
opts.restart = 10;

% Iterative refinement 
opts.nitref = 1;           % max # iterative refinement steps
opts.force_itref = false;  % force iterative refinement (true/false)
opts.itref_tol = 1.0e-12;
% opts.residual_update = true;

% diary('cpkrylov_unsym.txt')

% CPDQGMRES
fprintf('\n\n******************* %s mem %d *******************\n\n', cpk_string3, opts.mem)
[ cpk3x,  cpk3flags,  cpk3residHistory] =  reg_cpkrylov(cpk3, b, Q, A, C, G, opts);
cpk3x1  = cpk3x(1:n);
cpk3x2  = cpk3x(n+1:n+m);
fprintf('\n%s - rel err in (x, y) (refsol bsl): %7.1e', cpk_string3, norm(x-cpk3x)/norm(x));
fprintf('\n%s - rel err in x (refsol bsl): %7.1e', cpk_string3, norm(x(1:n)-cpk3x1)/norm(x(1:n)));
fprintf('\n%s - rel err in y (refsol bsl): %7.1e\n', cpk_string3, norm(x(n+1:n+m)-cpk3x2)/norm(x(n+1:n+m)));
fprintf('\n%s - iters = %d\n', cpk_string3, cpk3flags.niters);
pause

% restarted CPGMRES
fprintf('\n\n******************* %s(%d)*******************\n\n', cpk_string4, opts.restart)
[ cpk4x,  cpk4flags,  cpk4residHistory] =  reg_cpkrylov(cpk4, b, Q, A, C, G, opts);
cpk4x1  = cpk4x(1:n);
cpk4x2  = cpk4x(n+1:n+m);
fprintf('\n%s - rel err in (x, y) (refsol bsl): %7.1e', cpk_string4, norm(x-cpk4x)/norm(x));
fprintf('\n%s - rel err in x (refsol bsl): %7.1e', cpk_string4, norm(x(1:n)-cpk4x1)/norm(x(1:n)));
fprintf('\n%s - rel err in y (refsol bsl): %7.1e\n', cpk_string4, norm(x(n+1:n+m)-cpk4x2)/norm(x(n+1:n+m)));
fprintf('\n%s - iters = %d\n', cpk_string4, cpk4flags.niters);
pause

figure;
     
semilogy([0: size(cpk3residHistory.residHistory,1)-1], cpk3residHistory.residHistory, 'b.-', ...
         [0: size(cpk4residHistory.residHistory,1)-1], cpk4residHistory.residHistory, 'g:', ...
         'LineWidth', 2);

legend('CPDQGMRES', 'CPGMRES');
xlabel('iters');
ylabel('log(residual)');
title('Residual History');
