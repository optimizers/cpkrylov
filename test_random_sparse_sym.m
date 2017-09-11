% TEST PROGRAM
% Run constraint-preconditioned Krylov solvers on random symmetric
% generalized saddle-point systems.
%
% Daniela di Serafino, Sept. 10, 2017

% clear;
% close all;

n = 200; m = 100;
density = 0.01;
Q = sprandsym(n,density,1e-6,1);
% Q = sprand(n,n,density,1e-1);
% Q = Q + Q' + 1e5*speye(n);
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
% lambdaK = sort(real(eig(full(K))));
% n_lambdaK_pos = length(find(abs(lambdaK > 0)));
% n_lambdaK_neg = length(find(abs(lambdaK < 0)));
% fprintf('\n  num eigs > 0:  %d', n_lambdaK_pos);
% fprintf('\n  num eigs < 0:  %d\n', n_lambdaK_neg);
% plot(1:n+m,lambdaK,'*');
% pause

G = spdiags(abs(diag(Q)),0,n,n);
% G = speye(n);
P = [G  A' ; A  -C];
lambdaPK = sort(real(eig(full(K),full(P))));
n_lambdaPK_1 = length(find(abs(lambdaPK - 1) < 1e-8));
fprintf('\nCONSTRAINT-PRECONDITIONED MATRIX\n');
fprintf('\n  num. eigs close to 1: %d\n', n_lambdaPK_1);
figure
plot(1:n+m,lambdaPK,'*'); 
pause

cpk1 = @cpcg;
cpk_string1 = 'CPCG';
cpk2 = @cpcglanczos;
cpk_string2 = 'CPCGLANCZOS';
cpk3 = @cpdqgmres;
cpk_string3 = 'CPDQGMRES';
cpk4 = @cpgmres;
cpk_string4 = 'CPGMRES';

opts.print = true;     % true/false;
opts.atol = 1.0e-5;
opts.rtol = 1.0e-5;
opts.itmax = 600;
opts.mem = 3;          % 3 for symmetric systems;
opts.restart = 30;

% Iterative refinement 
opts.nitref = 1;          % max # iterative refinement steps
opts.force_itref = true;  % force iterative refinement (true/false)
% opts.itref_tol = 1.0e-15;
% opts.residual_update = true;

% diary('cpkrylov_sym.txt')

% CPCG
fprintf('\n\n******************* %s *******************\n\n', cpk_string1)
[ cpk1x,  cpk1flags,  cpk1residHistory] =  reg_cpkrylov(cpk1, b, Q, A, C, G, opts);
cpk1x1  = cpk1x(1:n);
cpk1x2  = cpk1x(n+1:n+m);
fprintf('\n%s - rel err in (x, y) (refsol bsl): %7.1e', cpk_string1, norm(x-cpk1x)/norm(x));
fprintf('\n%s - rel err in x (refsol bsl): %7.1e', cpk_string1, norm(x(1:n)-cpk1x1)/norm(x(1:n)));
fprintf('\n%s - rel err in y (refsol bsl): %7.1e\n', cpk_string1, norm(x(n+1:n+m)-cpk1x2)/norm(x(n+1:n+m)));
fprintf('\n%s - iters = %d\n', cpk_string1, cpk1flags.niters);
pause

% CPCGLANCZOS
fprintf('\n\n******************* %s *******************\n\n', cpk_string2)
[ cpk2x,  cpk2flags,  cpk2residHistory] =  reg_cpkrylov(cpk2, b, Q, A, C, G, opts);
cpk2x1  = cpk2x(1:n);
cpk2x2  = cpk2x(n+1:n+m);
fprintf('\n%s - rel err in (x, y) (refsol bsl): %7.1e', cpk_string2, norm(x-cpk2x)/norm(x));
fprintf('\n%s - rel err in x (refsol bsl): %7.1e', cpk_string2, norm(x(1:n)-cpk2x1)/norm(x(1:n)));
fprintf('\n%s - rel err in y (refsol bsl): %7.1e\n', cpk_string2, norm(x(n+1:n+m)-cpk2x2)/norm(x(n+1:n+m)));
fprintf('\n%s - iters = %d\n', cpk_string2, cpk2flags.niters);
pause

% CPDQGMRES
fprintf('\n\n******************* %s - mem %d *******************\n\n', cpk_string3, opts.mem)
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
      
semilogy([0: size(cpk1residHistory.residHistory,1)-1], cpk1residHistory.residHistory, 'm--', ...
         [0: size(cpk2residHistory.residHistory,1)-1], cpk2residHistory.residHistory, 'r-', ...
         [0: size(cpk3residHistory.residHistory,1)-1], cpk3residHistory.residHistory, 'b.-', ...
         [0: size(cpk4residHistory.residHistory,1)-1], cpk4residHistory.residHistory, 'g:', ...
         'LineWidth', 2);

legend('CPCG', 'CPCGLANCZOS', 'CPDQGMRES', 'CPGMRES');
xlabel('iters');
ylabel('log(residual)');
title('Residual History');
