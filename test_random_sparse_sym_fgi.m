% TEST PROGRAM
% Run constraint-preconditioned Krylov solvers on random symmetric
% generalized saddle-point systems.
%
% Daniela di Serafino, Sept. 10, 2017

% clear;
% close all;

n = 500; m = 200;

% Build Q.
densityQ = 0.01;
rcQ = 1.e-8;
Q = sprandsym(n,densityQ,rcQ,1);
nz = nnz(Q);
% Add alpha to nn diagonal entries of Q. Comment the following lines
% if you do not want to modify Q.
dq = 1e-12*ones(n,1);
alpha = 5.;
nn = n/4;
dq(1:floor(nn),1) = alpha;
Q = Q + spdiags(dq,0,n,n);
condest(Q)
% dd = spdiags(diag(Q),0,n,n);
% d1 = spdiags(ones(n,1),0,n,n);
% condest(dd)
% condest(Q-dd+d1)


% Build A.
rcA = 1.e-1;
densityA = densityQ;
A = sprand(m,n,densityA,rcA);
% A = 1e8*A;

% Build C.
rcC = 1.e-8;
density = 0.8;
beta = 1e-4;
% Different ways of building C. Uncomments your choice.
% C = spdiags(diag(sprandsym(m,0,rcC,1)),0,m,m);           % C diag spd
% C = spdiags(abs(diag(sprand(m,m,density,rcC))),0,m,m);   % C diag with entries >= 0
C = beta*speye(m);                                       % C = beta*I

K = [Q  A' ; A  -C];
b = rand(n+m, 1);

spy(K)

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
% set(gca, 'Fontsize', 14);
% title('Eigenvalues of M', 'Fontsize', 16);

% pause

G = spdiags(abs(diag(Q)),0,n,n);
% G = speye(n);
P = [G  A' ; A  -C];
lambdaPK = sort(real(eig(full(K),full(P))));
n_lambdaPK_1 = length(find(abs(lambdaPK - 1) < 1e-8));
fprintf('\nCONSTRAINT-PRECONDITIONED MATRIX\n');
fprintf('\n  num. eigs close to 1: %d\n', n_lambdaPK_1);
figure
semilogy(1:n+m,lambdaPK,'*');
set(gca, 'Fontsize', 14);
title('Eigenvalues of P^{-1}M', 'Fontsize', 16);
pause

cpk1 = @cpcg;
cpk_string1 = 'CPCG';
cpk22 = @cpcglanczos2;
% cpk_string22 = 'CPCGLANCZOS2';
% cpk2 = @cpcglanczos;
% cpk_string2 = 'CPCGLANCZOS';
cpk3 = @cpdqgmres;
cpk_string3 = 'CPDQGMRES';
cpk32 = @cpdqgmres;
cpk_string32 = 'CPDQGMRES';
% cpk4 = @cpgmres;
% cpk_string4 = 'CPGMRES';

opts.print = false;       % true/false;
opts.atol = 1.0e-6;
opts.rtol = 1.0e-6;
opts.itmax = 500;
opts.mem = 2;             % 2 for symmetric systems;
opts.restart = 20;

% Iterative refinement. 
opts.nitref = 1;          % max # iterative refinement steps
opts.force_itref = true;  % force iterative refinement (true/false)
% opts.itref_tol = 1.0e-15;
% opts.residual_update = true;

% diary('cpkrylov_sym.txt')

% CPCG.
fprintf('\n\n******************* %s *******************\n\n', cpk_string1)
[cpk1x,  cpk1stats,  cpk1flag] =  reg_cpkrylov(cpk1, b, Q, A, C, G, opts);
cpk1x1  = cpk1x(1:n);
cpk1x2  = cpk1x(n+1:n+m);
fprintf('\n%s - rel err in (x, y) (refsol bsl): %7.1e', cpk_string1, norm(x-cpk1x)/norm(x));
fprintf('\n%s - rel err in x (refsol bsl): %7.1e', cpk_string1, norm(x(1:n)-cpk1x1)/norm(x(1:n)));
fprintf('\n%s - rel err in y (refsol bsl): %7.1e\n', cpk_string1, norm(x(n+1:n+m)-cpk1x2)/norm(x(n+1:n+m)));
fprintf('\n%s - iters = %d\n', cpk_string1, cpk1stats.niters);
pause

% % CPCGLANCZOS2.
% fprintf('\n\n******************* %s *******************\n\n', cpk_string22)
% [cpk22x,  cpk22stats,  cpk22flag] =  reg_cpkrylov(cpk22, b, Q, A, C, G, opts);
% cpk22x1  = cpk22x(1:n);
% cpk22x2  = cpk22x(n+1:n+m);
% fprintf('\n%s - rel err in (x, y) (refsol bsl): %7.1e', cpk_string22, norm(x-cpk22x)/norm(x));
% fprintf('\n%s - rel err in x (refsol bsl): %7.1e', cpk_string22, norm(x(1:n)-cpk22x1)/norm(x(1:n)));
% fprintf('\n%s - rel err in y (refsol bsl): %7.1e\n', cpk_string22, norm(x(n+1:n+m)-cpk22x2)/norm(x(n+1:n+m)));
% fprintf('\n%s - iters = %d\n', cpk_string22, cpk22stats.niters);
% pause
% 
% % CPCGLANCZOS.
% fprintf('\n\n******************* %s *******************\n\n', cpk_string2)
% [cpk2x,  cpk2stats,  cpk2flag] =  reg_cpkrylov(cpk2, b, Q, A, C, G, opts);
% cpk2x1  = cpk2x(1:n);
% cpk2x2  = cpk2x(n+1:n+m);
% fprintf('\n%s - rel err in (x, y) (refsol bsl): %7.1e', cpk_string2, norm(x-cpk2x)/norm(x));
% fprintf('\n%s - rel err in x (refsol bsl): %7.1e', cpk_string2, norm(x(1:n)-cpk2x1)/norm(x(1:n)));
% fprintf('\n%s - rel err in y (refsol bsl): %7.1e\n', cpk_string2, norm(x(n+1:n+m)-cpk2x2)/norm(x(n+1:n+m)));
% fprintf('\n%s - iters = %d\n', cpk_string2, cpk2stats.niters);
% pause

% CPDQGMRES.
fprintf('\n\n******************* %s - mem %d *******************\n\n', cpk_string3, opts.mem)
[cpk3x,  cpk3stats,  cpk3flag] =  reg_cpkrylov(cpk3, b, Q, A, C, G, opts);
cpk3x1  = cpk3x(1:n);
cpk3x2  = cpk3x(n+1:n+m);
fprintf('\n%s - rel err in (x, y) (refsol bsl): %7.1e', cpk_string3, norm(x-cpk3x)/norm(x));
fprintf('\n%s - rel err in x (refsol bsl): %7.1e', cpk_string3, norm(x(1:n)-cpk3x1)/norm(x(1:n)));
fprintf('\n%s - rel err in y (refsol bsl): %7.1e\n', cpk_string3, norm(x(n+1:n+m)-cpk3x2)/norm(x(n+1:n+m)));
fprintf('\n%s - iters = %d\n', cpk_string3, cpk3stats.niters);
pause

opts.mem = 15; 

% CPDQGMRES.
fprintf('\n\n******************* %s - mem %d *******************\n\n', cpk_string32, opts.mem)
[cpk32x,  cpk32stats,  cpk32flag] =  reg_cpkrylov(cpk32, b, Q, A, C, G, opts);
cpk32x1  = cpk32x(1:n);
cpk32x2  = cpk32x(n+1:n+m);
fprintf('\n%s - rel err in (x, y) (refsol bsl): %7.1e', cpk_string32, norm(x-cpk32x)/norm(x));
fprintf('\n%s - rel err in x (refsol bsl): %7.1e', cpk_string32, norm(x(1:n)-cpk32x1)/norm(x(1:n)));
fprintf('\n%s - rel err in y (refsol bsl): %7.1e\n', cpk_string32, norm(x(n+1:n+m)-cpk32x2)/norm(x(n+1:n+m)));
fprintf('\n%s - iters = %d\n', cpk_string32, cpk32stats.niters);
pause

% % restarted CPGMRES.
% fprintf('\n\n******************* %s(%d)*******************\n\n', cpk_string4, opts.restart)
% [cpk4x,  cpk4stats,  cpk4flag] =  reg_cpkrylov(cpk4, b, Q, A, C, G, opts);
% cpk4x1  = cpk4x(1:n);
% cpk4x2  = cpk4x(n+1:n+m);
% fprintf('\n%s - rel err in (x, y) (refsol bsl): %7.1e', cpk_string4, norm(x-cpk4x)/norm(x));
% fprintf('\n%s - rel err in x (refsol bsl): %7.1e', cpk_string4, norm(x(1:n)-cpk4x1)/norm(x(1:n)));
% fprintf('\n%s - rel err in y (refsol bsl): %7.1e\n', cpk_string4, norm(x(n+1:n+m)-cpk4x2)/norm(x(n+1:n+m)));
% fprintf('\n%s - iters = %d\n', cpk_string4, cpk4stats.niters);
% pause

figure;
      
semilogy([0: size(cpk1stats.residHistory,1)-1],  cpk1stats.residHistory,  'r--', ...
         [0: size(cpk3stats.residHistory,1)-1],  cpk3stats.residHistory,  'b-.', ...
         'LineWidth', 2);
%          [0: size(cpk22stats.residHistory,1)-1], cpk22stats.residHistory, 'k:+', ...
%          [0: size(cpk2stats.residHistory,1)-1],  cpk2stats.residHistory,  'r-', ...
%          [0: size(cpk3stats.residHistory,1)-1],  cpk3stats.residHistory,  'b-.', ...
%          [0: size(cpk4stats.residHistory,1)-1],  cpk4stats.residHistory,  'g:', ...
%          'LineWidth', 2);

minx = 0;
maxx = max( size(cpk1stats.residHistory,1)-1, size(cpk3stats.residHistory,1)-1 );
miny = min( min(cpk1stats.residHistory), min(cpk3stats.residHistory) );
maxy = max( max(cpk1stats.residHistory), min(cpk3stats.residHistory) );
set(gca, 'Fontsize', 14);
axis([minx maxx miny maxy]);
xlabel('iters', 'FontSize', 16);
ylabel('log(||res||)', 'FontSize', 16);
title('residual norm history', 'FontSize', 16);
id = legend('CPCG', 'CPDQGMRES - mem = 2');
set(id, 'FontSize', 16);
%legend('CPCG', 'CPCGLANCZOS2', 'CPCGLANCZOS', 'CPDQGMRES', 'CPGMRES');

figure
semilogy([0: size(cpk1stats.residHistory,1)-1],  cpk1stats.residHistory,  'r--', ...
         [0: size(cpk3stats.residHistory,1)-1],  cpk3stats.residHistory,  'b-.', ...
         [0: size(cpk32stats.residHistory,1)-1], cpk32stats.residHistory, 'k-',  ...
         'LineWidth', 2);

minx = 0;
maxx = max( size(cpk1stats.residHistory,1)-1, size(cpk3stats.residHistory,1)-1 );
miny = min( min( min(cpk1stats.residHistory), min(cpk3stats.residHistory) ), min(cpk32stats.residHistory) );
maxy = max( max( max(cpk1stats.residHistory), max(cpk3stats.residHistory) ), max(cpk32stats.residHistory) );
set(gca, 'Fontsize', 14);
axis([minx maxx miny maxy]);
xlabel('iters', 'FontSize', 16);
ylabel('log(||res||)', 'FontSize', 16);
title('residual norm history', 'FontSize', 16);
id = legend('CPCG', 'CPDQGMRES - mem = 2', 'CPDQGMRES - mem = 15');
set(id, 'FontSize', 16);
