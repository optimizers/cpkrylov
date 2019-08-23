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
- use `logging4matlab`?
