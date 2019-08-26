# TO DO

- the iterative methods should not need B as input (they should always start from 0) --- DONE
- change driver so we always use 0 as starting point --- DONE
- move niters from flag to stats --- DONE
- main driver should just pass `b1`, no need to pass the block of zeros --- DONE
- the `dy` coming out of the iterative solver results in an incorrect `y` if we update `y = y0 + dy` --- FIXED
- CPMINRES --- DONE
- CPSYMMLQ --- DONE
- write a simple example program - DONE (2 example programs)
- write README for github
- which licence?

---------------------------------------------------------------------------
FUTURE WORK
- CP-GMRES, CP-DQGMRES: use Householder transformations in the generation
  of the Krylov basis, in order to reduce the numerical error?
  See H.F. Walker, Implementation of the GMRES Method Using Householder
  Transformations, SIAM J. Sci. Comp. 9(1), 1988.
- Add further stopping criteria? Backward error stopping criterion already
  implemented in CP-CGLANCZOS.
- Possible improvement of Lanczos process implementation.
- use `logging4matlab`?
