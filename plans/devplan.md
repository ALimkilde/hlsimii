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

## 2. Implement Static Mechanics kernel

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
- [X] Add pull_to_tension functionality

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
- [X] Compare RODAS4 against RK45 on representative problems.
- [X] Verify adaptive timestep behavior and convergence.

## 11. RODAS4 reference implementation

- [X] Make the dense RODAS4 implementation reliable.
- [X] Add regression tests.
- [X] Save representative trajectories/results.
- [X] Treat this as the reference implementation for later optimization.
- [X] Add event handling:
  - [X] Detect event crossings
  - [X] Locate roots
  - [X] Restart integration at events
- [X] Start with robust/simple root finding such as bisection before optimizing it.

## 12. Structured Jacobian solver

- [X] Exploit the mechanical Jacobian structure.
- [X] Use the state partition `y = (q, v)`.
- [X] Eliminate the position variables analytically from the Rosenbrock linear system.
- [X] Reduce the system to a 2N × 2N block-tridiagonal system.
- [X] Represent the reduced system using 2×2 blocks.
- [X] Implement block LU/Thomas factorization.
- [X] Implement a small, robust 2×2 solve rather than explicitly forming inverses.
- [X] Reuse all factorization/workspace memory.
- [X] Compare every structured solve against the dense reference solver.

## 13. Performance optimization

- [X] Profile before optimizing further.
- [X] Measure separately:
  - [X] Mechanics/RHS evaluation
  - [X] Jacobian construction
  - [X] Factorization
  - [X] Linear solves
  - [X] Event handling
  - [X] Total integration
- [-] Eliminate allocations from hot loops.
- [-] Reuse all temporary vectors and matrices.
- [X] Benchmark alternative small-matrix implementations.
- [X] Keep automated correctness tests alongside performance tests.
- [X] Profile for really large N and different lines

## 14. Found optimization

- [X] Assemble; reuse sqrt and unit vector and 1/len
- [-] Precompute; vlen?
- [X] Avoid rhsjac on rejected steps
- [X] Make event leash only compute two forces
- [X] Template edge_drag, edge_force -> 10% speedup on WASM
- [-] Avoid allocations in Rodas, physics -> 3% speedup
- [X] Precompute D[i] L[l] and store in L[l] -> 3% speedup
- [X] Rewrite drag derivative to rank-1

## 15. Make sure IO is correct

- [X] Make sure peak forces are calculated within each timestep (not only on output grid)
- [X] Make sure the whole table is computed correctly (implemented by me)
- [-] Save output to .csv file
- [X] Take Drag and KV scaling as input
- [X] Take Webbing width as input
- [-] Include Shockload values 
   - [-] Peak webbing during backup release
   - [-] Peak leash during flick up during backup release
- [X] Extensive testing against ISA backup fall

## 16. Energy calculation

- [X] Make 'static' energy solver for lower point of fall (assume f = 0 for all other nodes, fx = 0 for slackliner node, and Ep_old + Ee_old = Ep_new + Ee_new)
- [X] Test the solver
- [ ] Modify to shift in x-direction
- [ ] Make solver with mass lumped in the middle.
- [ ] Test against ISA - should give exactly the same!
- [ ] Reformulate the energy solver with the leash force P as the unknown (fixes the singular start: the energy row of the Jacobian is 0 at static equilibrium)
   - Inner: static solve with an extra downward force P on the slackliner's node (reuses static_solver and the block-tridiagonal LU)
   - Outer: solve the scalar equation E(q(P)) - E_goal = 0 with bisection or secant
   - P >= 0 picks the lower of the two solutions; P is the peak leash force
   - To figure out: how to pass the external force to the static solver, whether E(P) is monotonic, how to bracket P

## 16. WASM

- [X] Compile the existing correct implementation to WASM.
- [X] Keep the entire integration run inside WASM where possible.
- [X] Minimize JavaScript ↔ WASM calls.
- [X] Benchmark native C++ and WASM versions.
- [X] Enable compiler optimizations.
- [X] Investigate WASM SIMD only if profiling justifies it.
- [ ] Verify numerical results against the native/reference implementation.

## 17. Website integration

- [X] Define a small, stable WASM API.
- [X] Seperate backup and leash falls so they can be runned in parallel
- [X] Keep physics/integration state inside WASM.
- [X] Transfer only necessary inputs/results across the JS/WASM boundary.
- [X] Add visualization separately from the numerical core.
- [-] Add performance regression tests for representative simulations.
- [ ] Write test that web site agrees with hlsim

## 18. UI
- [X] Make run button follow header
- [X] Move run botton to the right
- [X] Make settings for slackliner fit more nicely
- [X] Make cookies based saves of setups
- [ ] Add a few examples
- [-] Reconfigure zoom
- [-] Make a not to scale button? (no axis equal)
- [X] Make documentation .pdf available
- [X] Make a warning about not using this instead of expert knowledge
- [ ] Rewrite all sentences
- [X] Make step 0.1 for most settings
- [X] Mechanism to reload/restart backend
- [X] Delete height curve
- [X] Move the kv and drag scale to bottom
- [X] Visualice harness and legs differently

## 19. Bugs

- [ ] Pull line cannot pull past a segment.
- [ ] Show impact (from leash taught frame max)? Or mark that it is not impact but max tension shown?

## 19. Validations

- [ ] Finish documentation!
- [ ] Validate and make report on experiment
- [ ] Validate and make report on comparison to ISA backup html

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
