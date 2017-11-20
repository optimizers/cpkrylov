% TEST PROGRAM
% Run CPDQGMRES and CPGMRES(r) on random unsymmetric regularized
% saddle-point systems.
%
% Daniela di Serafino, Nov. 20, 2017

% clear;
% close all;

n = 300; m = 100;
density = 0.01;
rc = 1e-8;
Q = sprand(n,n,density,rc);    % Q may have zeros on diag
a = 5.0e-1;                    % the smaller a, the worse Q.
Q = Q + a*speye(n);
A = sprand(m,n,density,rc);
% A = 1e8*A;
% C = spdiags(diag(sprandsym(m,0,rc,1)),0,m,m);           % C diag spd
density2 = 0.8;
C = spdiags(abs(diag(sprand(m,m,density2,rc))),0,m,m);  % C diag with entries >= 0
% C = 1e-2*C;
K = [Q  A' ; A  -C];
b = rand(n+m, 1);

x = K\b;

fprintf('\nSADDLE-POINT SYSTEM (n = %d, m = %d)\n', n, m);
fprintf('\n  [ A   B'' ] [x] = [b]');
fprintf('\n  [ B  -C  ] [y]   [0]\n');
fprintf('\n  matrix rank: %d\n', rank(full(K)));

% Info on the eigenvalues of the saddle-point matrix.
lambdaK = eig(full(K));
figure
plot(lambdaK,'*');
title('Eigenvalues of saddle-point matrix');

pause

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
fprintf('\n  num. eigs Inf (there should not be any): %d\n', n_lambdaPK_inf);
figure
plot(lambdaPK,'*');
title('Eigenvalues of preconditioned matrix');

pause

cpk3 = @cpdqgmres;
cpk_string3 = 'CPDQGMRES';
cpk4 = @cpgmres;
cpk_string4 = 'CPGMRES';
% cpk5 = @cpdqgmres_old;
% cpk_string5 = 'CPDQGMRES_OLD';

opts.print = false;     % true/false;
opts.atol = 1.0e-6;
opts.rtol = 1.0e-6;
opts.itmax = 500;
opts.mem = 10;
opts.restart = 10;

% Iterative refinement 
opts.nitref = 1;           % max # iterative refinement steps
opts.force_itref = true;  % force iterative refinement (true/false)
opts.itref_tol = 1.0e-12;
% opts.residual_update = true;

% diary('cpkrylov_unsym.txt')

% CPDQGMRES
fprintf('\n\n******************* %s mem %d *******************\n\n', cpk_string3, opts.mem)
[ cpk3x,  cpk3stats, cpk3flag] =  reg_cpkrylov(cpk3, b, Q, A, C, G, opts);
cpk3x1  = cpk3x(1:n);
cpk3x2  = cpk3x(n+1:n+m);
fprintf('\n%s - rel err in (x, y) (refsol bsl): %7.1e', cpk_string3, norm(x-cpk3x)/norm(x));
fprintf('\n%s - rel err in x (refsol bsl): %7.1e', cpk_string3, norm(x(1:n)-cpk3x1)/norm(x(1:n)));
fprintf('\n%s - rel err in y (refsol bsl): %7.1e\n', cpk_string3, norm(x(n+1:n+m)-cpk3x2)/norm(x(n+1:n+m)));
fprintf('\n%s - iters = %d\n', cpk_string3, cpk3stats.niters);
pause

% restarted CPGMRES
fprintf('\n\n******************* %s(%d)*******************\n\n', cpk_string4, opts.restart)
[ cpk4x,  cpk4stats, cpk4flag] =  reg_cpkrylov(cpk4, b, Q, A, C, G, opts);
cpk4x1  = cpk4x(1:n);
cpk4x2  = cpk4x(n+1:n+m);
fprintf('\n%s - rel err in (x, y) (refsol bsl): %7.1e', cpk_string4, norm(x-cpk4x)/norm(x));
fprintf('\n%s - rel err in x (refsol bsl): %7.1e', cpk_string4, norm(x(1:n)-cpk4x1)/norm(x(1:n)));
fprintf('\n%s - rel err in y (refsol bsl): %7.1e\n', cpk_string4, norm(x(n+1:n+m)-cpk4x2)/norm(x(n+1:n+m)));
fprintf('\n%s - iters = %d\n', cpk_string4, cpk4stats.niters);
pause

% opts.mem = 11;
% 
% % CPDQGMRES_OLD
% fprintf('\n\n******************* %s mem %d *******************\n\n', cpk_string5, opts.mem)
% [ cpk5x,  cpk5stats,  cpk5flag] =  reg_cpkrylov(cpk5, b, Q, A, C, G, opts);
% cpk5x1  = cpk5x(1:n);
% cpk5x2  = cpk5x(n+1:n+m);
% fprintf('\n%s - rel err in (x, y) (refsol bsl): %7.1e', cpk_string5, norm(x-cpk5x)/norm(x));
% fprintf('\n%s - rel err in x (refsol bsl): %7.1e', cpk_string5, norm(x(1:n)-cpk5x1)/norm(x(1:n)));
% fprintf('\n%s - rel err in y (refsol bsl): %7.1e\n', cpk_string5, norm(x(n+1:n+m)-cpk5x2)/norm(x(n+1:n+m)));
% fprintf('\n%s - iters = %d\n', cpk_string5, cpk5stats.niters);
% pause

figure;
     
semilogy([0: size(cpk3stats.residHistory,1)-1], cpk3stats.residHistory, 'b.-', ...
         [0: size(cpk4stats.residHistory,1)-1], cpk4stats.residHistory, 'g:', ...
         'LineWidth', 2);
%          [0: size(cpk5stats.residHistory,1)-1], cpk5stats.residHistory, 'r--', ...
%          'LineWidth', 2);

legend('CPDQGMRES', 'CPGMRES');
% legend('CPDQGMRES', 'CPGMRES', 'CPDQGMRES-OLD');
xlabel('iters');
ylabel('log(residual)');
title('Residual History');
