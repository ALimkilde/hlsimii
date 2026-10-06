# Structure decisions

- **Slackliner (static):** a `PointMass{node, m}` passed as an argument (`= nullptr` by default). It's never stored in `LineModel`.
- **Placement:** `nearest_node_x(q, x)` is a free function. The experiment does: empty solve → nearest node → loaded solve, warm-started from `q0` (with load stepping if needed).
- **Dynamics:** a hybrid system with a rigid leash. Free fall is analytic, the impact is an inelastic map, and when attached the node carries `m_k + m_s`. `y` holds the line only, so the Jacobian stays block-tridiagonal.
- **State vs. configuration:** state lives in `y`. Phase and configuration reach `rhs` through a lambda, so `LineModel` stays `const`.
- **Layers:** experiment → hybrid driver → solvers (Newton, Rodas4) → `LineModel` → mesh → line. Each layer knows only the ones below it.
- **pull_to_tension:** lives in the experiment layer, built from `anchor_tension` (physics), a generic `find_root` and a lambda. Rebuilding `DiscreteLine` every iteration is fine.
- **Main/backup:** keep one model with shared nodes. Switch to a loop over the two springs when you add damping.
- **Mainline break:** `with_broken_mainline(DiscreteLine)` returns a copy with `k_main = 0` on every element, and the experiment builds a second `LineModel` from it. Watch `any_edge_slack`: it should test for zero tension, not lengths.
- **static_solver:** move it out of `LineModel` as a free function, with its settings (`maxsteps`, `reduce_alpha`, the tolerance) as options. The initial guess stays in `LineModel`.
- **mesh.h / physics.h:** keep them separate. Write down the node-index mapping once in `physics.h`.
- **Warm start for pull_to_tension:** not needed.
