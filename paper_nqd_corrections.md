# NQD paper-content corrections (from implementation validation)

This note records the exact issue we found in the current write-up and the precise correction.

## What was wrong

In the current text, the closed-form equations for
`sum_{u in U_xy} c(u,x)` and `sum_{u in U_xy} c(u,y)` include:

- `D(x) - D(parent(w))` and `D(y) - D(parent(w))`, and
- `P(y) - P(w)` and `P(x) - P(w)`.

Where `w = LCA(x,y)`.

This is not consistent with the actual node set `U_xy` (internal nodes on the x-y path excluding endpoints), because:

1. Using `D(parent(w))` includes the `w` contribution when it should be excluded from `U_xy`.
2. Using `P(y) - P(w)` includes the leaf node `y` term `p(y)` (and symmetrically includes `p(x)`), but `y` and `x` are endpoints and are excluded from `U_xy`.

These two inclusions caused mismatches against the baseline path-traversal implementation.

## Correct formulas

Let `par(z)` be parent of `z` and treat prefix values at an undefined parent as 0.

The corrected sums are:

- `sum_{u in U_xy} c(u,x) = ( D(x) - D(w) ) + ( P(par(y)) - P(w) )`
- `sum_{u in U_xy} c(u,y) = ( D(y) - D(w) ) + ( P(par(x)) - P(w) )`

And still:

- `|U_xy| = depth(x) + depth(y) - 2*depth(w) - 1`
- `NQD_gt(x,y) = 1/2 * ( sum_cx + sum_cy - 2*|U_xy| )`

## Edge-case conventions

- If `x = y`, return 0.
- If `par(x)` or `par(y)` is undefined (degenerate cases), treat corresponding prefix as 0.
- If `|U_xy| < 0` (tiny-tree degeneracy), clamp to 0.

## Validation status

After applying the corrected formulas, optimized and baseline methods match in validation mode (worst diff = 0 on tested inputs), including final species distance matrix comparison.
