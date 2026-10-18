# Rodas4 implementation plan

- [X] Implement the coefficient and stages
- [X] Decide on tiny test problem - CLAUDE did this
- [X] Get a fixed stepsize version running
- [ ] Add statistics counters
- [X] Add stepsize control and error estimation
- [X] Add dense output and observers
- [X] Add event handling
- [X] Make linalg solver abstract
- [X] Implement the full gradients of physics rhs
- [X] Implement framework for getting jacobian of rhs.
- [ ] Test on this simulation repo's rhs
- [ ] Make sure that jacobian is not regenerated after rejected steps
- [ ] Make sure that we don't have unecessary rhs calls
- [ ] Make sure that we don't have unecessary memory (ynew? embedded?)
- [ ] Exploit structures of Jacobian
- [ ] Try the predictive stepsize constoller from Gustafson (see fortran code)
