# Webbing stiffness: from secant to stretch curves

## Problem
Stiffness today comes from one point: `k = 100·T / (L·stretch_pct)`. That is the
*secant* stiffness from zero, which includes the soft "toe" where the weave
straightens. A highline never goes near zero tension. It works from pretension
upward, so the slope that matters is the one above pretension (the *tangent*).
Using the secant makes lines too soft. In the ISA 2024 validation, nylon needed
`--stiffness-scale ≈ 1.5`.

## Three effects behind the gap
1. **Secant vs. tangent.** Geometry of the stretch curve.
2. **Load history.** A line that has been rigged a while has bedded in and crept,
   and on reloading it follows a steeper curve than on first loading. This is
   the known "retension after hanging" and "permanent lines lose their bounce,
   rest the webbing to recover" effect.
3. **Rate dependence.** Nylon and polyester are viscoelastic and respond stiffer
   to fast loading. Probably a modest effect (roughly logarithmic in rate).

We can't estimate 2 and 3 for a given user's line, so we don't try.

## Approach
- **Pretension (or sag) is a user input.** It is measured, not derived from slack
  length, so creep's static effect (longer line, lower tension) is absorbed.
  *To check: what the sim currently takes as input.*
- **User enters the data-sheet stretch curve**: 1–N points (kN → %), ideally
  covering pretension up to typical peaks (~2–10 kN). Plus the material class
  (nylon / polyester / Dyneema). One point is a backward-compatible fallback.
- **The simulator brackets the stiffness:**
  - *Fresh / rested*: the data-sheet curve as given (soft end).
  - *Settled*: the curve's slope × a per-material settled factor (stiff end),
    e.g. nylon ~1.5, Dyneema ~1. Calibrated by us, never entered by users.
- **Report outputs as a range**, since each end is the worst case for something
  different:

  | Quantity                     | Worst case |
  |------------------------------|------------|
  | Peak anchor / leash forces   | settled (stiff) |
  | Fall clearance / lowest point| fresh (soft) |

## Display units
The physical quantity is EA in kN, which means nothing to users. Show the
derived **%/kN at pretension** ("each extra kN adds x %") as a read-only value
next to the curve.

## Catalogue row
`name | material | stretch points (kN → %) | g/m | width`

## Open questions
- Nonlinear curve in the solver, or linearize around pretension? The second is a
  smaller change.
- Calibrating settled factors: bounce frequency on long-rigged lines vs. the
  data-sheet slope; unload/reload test curves; check that one factor fits all ISA
  trials (if the best factor grows with fall severity, rate matters).
