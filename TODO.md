# TO DO

- the iterative methods should not need B as input (they should always start from 0) --- DONE
- change driver so we always use 0 as starting point --- DONE
- use `logging4matlab`
- main driver should just pass `b1`, no need to pass the block of zeros --- DONE
- the `dy` coming out of the iterative solver results in an incorrect `y` if we update `y = y0 + dy` --- FIXED for all solvers but cpcg
- MINRES
