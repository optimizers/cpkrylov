# TO DO (May 30, 2020)

- Add further stopping criteria? Backward error stopping criterion already
  implemented in CP-CGLANCZOS.
- Possible improvement of Lanczos process implementation, for better stability.
- CP-GMRES, CP-DQGMRES: use Householder transformations in the generation
  of the Krylov basis, in order to reduce the numerical error.
  See H.F. Walker, Implementation of the GMRES Method Using Householder
  Transformations, SIAM J. Sci. Comp. 9(1), 1988.
- CPDQGMRES: reduce memory used for H.
- Write simpler example programs?
- Use `logging4matlab`?
