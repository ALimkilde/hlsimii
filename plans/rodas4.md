# Rodas4 implementation plan

## Implementations

- [X] Implement the coefficient and stages
- [X] Decide on tiny test problem - CLAUDE did this
- [X] Get a fixed stepsize version running
- [X] Add statistics counters
- [X] Add stepsize control and error estimation
- [X] Add dense output and observers
- [X] Add event handling
- [X] Make linalg solver abstract
- [X] Implement the full gradients of physics rhs
- [X] Implement framework for getting jacobian of rhs.
- [X] Test on this simulation repo's rhs

## Optimizations

- [ ] Make sure that jacobian is not regenerated after rejected steps
- [ ] Avoid unecessary rhs calls
- [ ] Avoid unecessary memory (ynew? embedded?)
- [ ] Avoid unecessary heap memory (stages to tmp_rhs)
- [ ] Exploit structures of Jacobian
   - [X] Exploit structures q | v structure => 5x speedup
   - [ ] Make subtracting (-operator) for BlockTriMat
   - [ ] Make MV product for BlockTriMat
   - [ ] Make add to diagonal method
- [ ] Try the predictive stepsize constoller from Gustafson (see fortran code)
