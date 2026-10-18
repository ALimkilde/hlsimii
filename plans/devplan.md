# Rope ODE Solver — Development Workflow

## 1. State representation

- [X] Define the highline segmented representation
- [X] Decide and document boundary representation.
- [X] Define the ODE state and indexing.
- [X] Use contiguous storage:
  - `q = [x0, y0, x1, y1, ...]`
  - `v = [vx0, vy0, vx1, vy1, ...]`
  - `state = [q | v]`
- [X] Use `Eigen::VectorXd` for the underlying state.
- [X] Use `Eigen::Map` for zero-copy N × 2 views of positions and velocities.

## 2. Implement Statix Mechanics kernel

- [X] Implement spring forces/accelerations with explicit loops.
- [X] Avoid unnecessary temporary allocations.
- [X] Keep mechanics independent of the ODE integrator.
- [X] Implement boundary conditions cleanly.
- [X] Add basic physics tests:
  - [X] Single spring
  - [X] Equilibrium
  - [X] Equal/opposite forces
  - [X] Gravity
  - [X] Slack/taut spring
  - [X] Boundary behavior

## 3. Static solver

- [X] Implement static equilibrium solves.
- [X] Use the same force calculation as the dynamic model.
- [X] Verify solutions physically.
- [X] Add slackliner to static equations
- [ ] Add pull_to_tension functionality <------ TODO

## 4. Spring Jacobian

- [X] Derive and implement the analytic spring Jacobian.
- [X] Use it to accelerate Newton/static solves.
- [X] Compare analytic derivatives against finite differences.
- [X] Test the Jacobian at many random physically reasonable configurations.
- [X] Investigate carefully around `max(..., 0)` transitions.

## 5. Extend Static Mechanics to Dynamic;

- [X] Refactor assemble to loop over edges. 
- [X] Make stuff from assemble reusable?
- [X] Implement rhs
- [X] Implement class for full fall simulations - header
- [X] Test on leash fall
- [X] Extend to backup fall
- [X] Figure out how to log data from the simulations!

## 5. RK45 reference integrator

- [X] Have Claude implement a clean adaptive RK45.
- [X] Keep the integrator independent of the mechanics.
- [X] Use RK45 as a reference solver rather than the final performance target.

## 6. Springs + RK45

- [X] Run the spring-only dynamics.
- [X] Compare against known/simple cases.

## 7. Add drag

- [X] Add drag to the mechanics kernel.
- [X] Validate against simple cases with known behavior.
- [X] Re-run RK45 tests.

## 8. Add Kelvin–Voigt damping

- [X] Add Kelvin–Voigt damping.
- [X] Validate independently.
- [X] Re-run the mechanical tests.
- [X] Compare behavior with/without damping.

## 9. Full analytic Jacobian

- [X] Derive the Jacobian of the complete system:
  - [X] Springs
  - [X] Drag
  - [X] Kelvin–Voigt damping
- [X] Validate against finite differences.
- [X] Keep the Jacobian implementation independent from the eventual optimized solver.

## 10. RODAS4 — first implementation

- [X] Verify the RODAS4 coefficients and equations.
- [X] Figure out a good interface for rhs and linear problems
- [X] Implement RODAS4 incrementally with Claude as a coding/review assistant.
  - [X] Fixed step size
  - [X] Adaptive step size
- [X] Initially use:
  - [X] Dense Jacobian
  - [X] Eigen dense LU
  - [ ] No structural optimization
- [ ] Compare RODAS4 against RK45 on representative problems.
- [ ] Verify adaptive timestep behavior and convergence.

## 11. RODAS4 reference implementation

- [ ] Make the dense RODAS4 implementation reliable.
- [ ] Add regression tests.
- [ ] Save representative trajectories/results.
- [ ] Treat this as the reference implementation for later optimization.
- [X] Add event handling:
  - [X] Detect event crossings
  - [X] Locate roots
  - [X] Restart integration at events
- [X] Start with robust/simple root finding such as bisection before optimizing it.

## 12. Structured Jacobian solver

- [ ] Exploit the mechanical Jacobian structure.
- [ ] Use the state partition `y = (q, v)`.
- [ ] Eliminate the position variables analytically from the Rosenbrock linear system.
- [ ] Reduce the system to a 2N × 2N block-tridiagonal system.
- [ ] Represent the reduced system using 2×2 blocks.
- [ ] Implement block LU/Thomas factorization.
- [ ] Implement a small, robust 2×2 solve rather than explicitly forming inverses.
- [ ] Reuse all factorization/workspace memory.
- [ ] Compare every structured solve against the dense reference solver.

## 13. Performance optimization

- [ ] Profile before optimizing further.
- [ ] Measure separately:
  - [ ] Mechanics/RHS evaluation
  - [ ] Jacobian construction
  - [ ] Factorization
  - [ ] Linear solves
  - [ ] Event handling
  - [ ] Total integration
- [ ] Eliminate allocations from hot loops.
- [ ] Reuse all temporary vectors and matrices.
- [ ] Benchmark alternative small-matrix implementations.
- [ ] Keep automated correctness tests alongside performance tests.

## 14. WASM

- [ ] Compile the existing correct implementation to WASM.
- [ ] Keep the entire integration run inside WASM where possible.
- [ ] Minimize JavaScript ↔ WASM calls.
- [ ] Benchmark native C++ and WASM versions.
- [ ] Enable compiler optimizations.
- [ ] Investigate WASM SIMD only if profiling justifies it.
- [ ] Verify numerical results against the native/reference implementation.

## 15. Website integration

- [ ] Define a small, stable WASM API.
- [ ] Keep physics/integration state inside WASM.
- [ ] Transfer only necessary inputs/results across the JS/WASM boundary.
- [ ] Add visualization separately from the numerical core.
- [ ] Add performance regression tests for representative simulations.

## 16. Energy calculation

- [ ] Make 'static' energy solver for lower point of fall (assume f = 0 for all other nodes, fx = 0 for slackliner node, and Ep_old + Ee_old = Ep_new + Ee_new)
- [ ] Test the solver
- [ ] Make solver with mass lumped in the middle.
- [ ] Test against ISA - should give exactly the same!

## Milestones

**Milestone 1 — Physics:**
Steps 1–4. Static equilibrium and analytic Jacobian are trustworthy.

**Milestone 2 — Reference solver:**
Steps 5–9. RK45 and the complete mechanics/Jacobian are validated.

**Milestone 3 — Working RODAS4:**
Steps 10–11. Dense RODAS4 is correct and serves as the reference implementation.

**Milestone 4 — Fast solver:**
Step 12. Structured 2×2 block-tridiagonal solver matches the dense implementation.

**Milestone 5 — Production:**
Steps 13–15. Profiled, optimized, compiled to WASM, and integrated into the website.
