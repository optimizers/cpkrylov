%     K2 - CPCG vs CPMINRES (atol = rtol = 1.00e-06, maxit = 2000, force_itref = 1)
% 
%         CPCG                                    CPMINRES
% 
% Problem   Iters     Tsetup     Tsolve       Ttot  Iters     Tsetup     Tsolve       Ttot
% 
%                                      IPiter = 0
%
dati_K2 = [
% cvxqp1_s
     62   6.52e-02   5.84e-02   1.36e-01     60   3.10e-02   2.90e-02   6.21e-02
% cvxqp1_m
    198   3.97e-02   2.04e-01   2.46e-01    176   2.04e-02   2.06e-01   2.30e-01
% cvxqp1_l
    546   7.29e-02   2.93e+00   3.01e+00    425   6.27e-02   2.29e+00   2.36e+00
% cvxqp2_s
     67   6.49e-03   4.56e-02   5.49e-02     66   9.69e-03   4.44e-02   5.48e-02
% cvxqp2_m
    222   1.10e-02   2.03e-01   2.16e-01    198   1.31e-02   1.76e-01   1.92e-01
% cvxqp2_l
    645   4.26e-02   3.00e+00   3.04e+00    525   4.30e-02   2.52e+00   2.56e+00
% cvxqp3_s
     52   1.06e-02   4.20e-02   5.63e-02     50   7.38e-03   2.85e-02   3.82e-02
% cvxqp3_m
    168   1.85e-02   1.55e-01   1.75e-01    142   1.44e-02   1.32e-01   1.48e-01
% cvxqp3_l
    471   8.33e-02   2.68e+00   2.77e+00    353   8.12e-02   2.06e+00   2.15e+00
% mosarqp1
      4   1.67e-02   2.17e-02   4.03e-02      4   1.69e-02   2.27e-02   4.19e-02
% mosarqp2
      7   1.09e-02   1.81e-02   3.15e-02      7   1.22e-02   2.14e-02   3.55e-02
%   stcqp1
     29   3.96e-02   9.42e-02   1.37e-01     28   2.97e-02   7.52e-02   1.06e-01
%   stcqp2
     37   2.68e-02   1.00e-01   1.29e-01     37   2.56e-02   9.50e-02   1.22e-01
%
%                                          IP it 5
%
% cvxqp1_s
     29   7.66e-02   5.86e-02   1.48e-01     29   3.78e-02   2.67e-02   6.72e-02
% cvxqp1_m
    160   3.88e-02   1.03e+00   1.07e+00    150   1.91e-02   9.78e-01   9.99e-01
% cvxqp1_l
    860   2.05e-01   5.29e+01   5.32e+01    737   2.14e-01   4.53e+01   4.55e+01
% cvxqp2_s
     37   4.41e-03   5.96e-02   6.64e-02     37   4.11e-03   4.63e-02   5.23e-02
% cvxqp2_m
    295   7.93e-03   1.65e+00   1.66e+00    273   6.98e-03   1.52e+00   1.53e+00
% cvxqp2_l
   1287   4.27e-02   6.95e+01   6.95e+01   1053   4.80e-02   5.74e+01   5.74e+01
% cvxqp3_s
     18   9.16e-03   3.38e-02   4.67e-02     18   4.58e-03   3.55e-02   4.28e-02
% cvxqp3_m
    109   2.26e-02   7.30e-01   7.56e-01    106   1.98e-02   6.81e-01   7.02e-01
% cvxqp3_l
    543   5.34e-01   3.49e+01   3.54e+01    416   5.10e-01   2.79e+01   2.84e+01
% mosarqp1
      3   2.02e-02   3.26e-02   5.67e-02      3   2.36e-02   2.73e-02   5.28e-02
% mosarqp2
      7   1.42e-02   7.66e-02   9.32e-02      7   1.10e-02   5.40e-02   6.73e-02
%   stcqp1
     24   2.31e-01   7.02e-01   9.37e-01     24   2.21e-01   7.50e-01   9.75e-01
%   stcqp2
     26   3.16e-02   7.35e-01   7.69e-01     26   3.00e-02   7.00e-01   7.32e-01
%
%                                         IP it 10
%
% cvxqp1_s
      9   7.42e-02   4.35e-02   1.31e-01      9   3.24e-02   2.15e-02   5.59e-02
% cvxqp1_m
     56   3.64e-02   4.71e-01   5.11e-01     54   1.83e-02   3.47e-01   3.67e-01
% cvxqp1_l
    289   2.16e-01   1.71e+01   1.73e+01    258   2.00e-01   1.54e+01   1.56e+01
% cvxqp2_s
      9   6.24e-03   1.27e-02   2.20e-02      9   7.85e-03   1.13e-02   2.12e-02
% cvxqp2_m
     84   8.98e-03   4.34e-01   4.45e-01     80   1.13e-02   4.18e-01   4.31e-01
% cvxqp2_l
    664   4.48e-02   3.60e+01   3.60e+01    569   4.71e-02   3.12e+01   3.12e+01
% cvxqp3_s
      6   9.30e-03   7.51e-03   1.96e-02      6   6.95e-03   6.96e-03   1.48e-02
% cvxqp3_m
     28   2.71e-02   1.65e-01   1.94e-01     28   2.04e-02   1.84e-01   2.06e-01
% cvxqp3_l
    129   5.31e-01   8.10e+00   8.63e+00    119   5.06e-01   7.62e+00   8.13e+00
%   stcqp1
     16   2.33e-01   5.06e-01   7.43e-01     16   2.21e-01   5.14e-01   7.38e-01 ];

iters_cg      = dati_K2(1:end,1);
tsetup_cg     = dati_K2(1:end,2);
tsolve_cg     = dati_K2(1:end,3);
ttot_cg       = dati_K2(1:end,4);
iters_minres  = dati_K2(1:end,5);
tsetup_minres = dati_K2(1:end,6);
tsolve_minres = dati_K2(1:end,7);
ttot_minres   = dati_K2(1:end,8);

iters=[iters_cg iters_minres];
figure;
perf_prof(iters,0)
% set(gca,'XTick',1:2:10)
%set(gca,'YTick',0:0.1:1)
xlabel('\chi', 'FontSize', 15);
ylabel('\pi(\chi)', 'FontSize', 15);
title('Iterations', 'FontSize', 14); 
legend('CP-CG', 'CP-MINRES', 'Location', 'SouthEast', 'FontSize', 14);
filefig = 'cg_vs_minres_iters.eps';
print(filefig, '-depsc');
close;

ttot=[ttot_cg ttot_minres];
figure;
perf_prof(ttot,0)
% set(gca,'XTick',1:2:10)
%set(gca,'YTick',0:0.1:1)
xlabel('\chi', 'FontSize', 15) ;
ylabel('\pi(\chi)', 'FontSize', 15);
title('Total time', 'FontSize', 14); 
legend('CP-CG', 'CP-MINRES', 'Location', 'SouthEast', 'FontSize', 14);
close;

tsolve=[tsolve_cg tsolve_minres];
figure;
perf_prof(tsolve,0)
% set(gca,'XTick',1:2:10)
%set(gca,'YTick',0:0.1:1)
xlabel('\chi', 'FontSize', 15);
ylabel('\pi(\chi)', 'FontSize', 15);
title('Solve time', 'FontSize', 14);
legend('CP-CG', 'CP-MINRES', 'Location', 'SouthEast', 'FontSize', 14);
filefig = 'cg_vs_minres_solve_time.eps';
print(filefig, '-depsc');
close;
