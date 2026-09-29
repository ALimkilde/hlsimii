# Rope ODE Solver — Development Workflow

## 1. State representation

- [ ] Define the ODE state and indexing.
- [ ] Use contiguous storage:
  - `q = [x0, y0, x1, y1, ...]`
  - `v = [vx0, vy0, vx1, vy1, ...]`
  - `state = [q | v]`
- [ ] Use `Eigen::VectorXd` for the underlying state.
- [ ] Use `Eigen::Map` for zero-copy N × 2 views of positions and velocities.
- [ ] Decide and document boundary representation.
- [ ] Define the highline segmented representation

## 2. Mechanics kernel

- [ ] Implement spring forces/accelerations with explicit loops.
- [ ] Avoid unnecessary temporary allocations.
- [ ] Keep mechanics independent of the ODE integrator.
- [ ] Implement boundary conditions cleanly.
- [ ] Add basic physics tests:
  - [ ] Single spring
  - [ ] Equilibrium
  - [ ] Equal/opposite forces
  - [ ] Gravity
  - [ ] Slack/taut spring
  - [ ] Boundary behavior

## 3. Static solver

- [ ] Implement static equilibrium solves.
- [ ] Use the same force calculation as the dynamic model.
- [ ] Verify solutions physically.
- [ ] Test difficult configurations and nearly slack springs.

## 4. Spring Jacobian

- [ ] Derive and implement the analytic spring Jacobian.
- [ ] Use it to accelerate Newton/static solves.
- [ ] Compare analytic derivatives against finite differences.
- [ ] Test the Jacobian at many random physically reasonable configurations.
- [ ] Investigate carefully around `max(..., 0)` transitions.

## 5. RK45 reference integrator

- [ ] Have Claude implement a clean adaptive RK45.
- [ ] Keep the integrator independent of the mechanics.
- [ ] Establish tolerances, error control, and output conventions.
- [ ] Use RK45 as a reference solver rather than the final performance target.

## 6. Springs + RK45

- [ ] Run the spring-only dynamics.
- [ ] Compare against known/simple cases.
- [ ] Check energy behavior where appropriate.
- [ ] Establish baseline performance and accuracy.

## 7. Add drag

- [ ] Add drag to the mechanics kernel.
- [ ] Validate against simple cases with known behavior.
- [ ] Re-run RK45 tests.

## 8. Add Kelvin–Voigt damping

- [ ] Add Kelvin–Voigt damping.
- [ ] Validate independently.
- [ ] Re-run the mechanical tests.
- [ ] Compare behavior with/without damping.

## 9. Full analytic Jacobian

- [ ] Derive the Jacobian of the complete system:
  - [ ] Springs
  - [ ] Drag
  - [ ] Kelvin–Voigt damping
- [ ] Validate against finite differences.
- [ ] Pay particular attention to signs and mass scaling.
- [ ] Keep the Jacobian implementation independent from the eventual optimized solver.

## 10. RODAS4 — first implementation

- [ ] Verify the RODAS4 coefficients and equations.
- [ ] Implement RODAS4 incrementally with Claude as a coding/review assistant.
- [ ] Initially use:
  - [ ] Dense Jacobian
  - [ ] Eigen dense LU
  - [ ] No structural optimization
- [ ] Compare RODAS4 against RK45 on representative problems.
- [ ] Verify adaptive timestep behavior and convergence.

## 11. RODAS4 reference implementation

- [ ] Make the dense RODAS4 implementation reliable.
- [ ] Add regression tests.
- [ ] Save representative trajectories/results.
- [ ] Treat this as the reference implementation for later optimization.
- [ ] Add event handling:
  - [ ] Detect event crossings
  - [ ] Locate roots
  - [ ] Restart integration at events
- [ ] Start with robust/simple root finding such as bisection before optimizing it.

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
