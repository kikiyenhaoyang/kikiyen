# GRIN Lie–Poisson regression tests

The standalone production page is Web/poisson_variety_grin.html. There are no
runtime build dependencies. Tests read and execute that page's actual inline
JavaScript; the numerical model is not copied into the test suite.

## Physical domain

For the UI's domain α>0, H>0, n₀>0, C≥0, let b=n₀²−H². The physical Gram matrix
requires X≥0, Y≥0, XY−Z²=C. Energy gives Y=b−αX, hence Z²=F(X)=−αX²+bX−C.

| Conditions | Physical energy–Casimir intersection |
| --- | --- |
| b<0 | Empty, even if the discriminant is nonnegative |
| b=0, C>0 | Empty |
| b=0, C=0 | Stationary vertex (0,0,0) |
| b>0, C=0 | Singular cone orbit, 0≤X≤b/α |
| b>0, C>0, Δ>0 | Guided orbit, X₋≤X≤X₊ |
| b>0, C>0, Δ=0 | Single reduced equilibrium (b/(2α),b/2,0) |
| b>0, C>0, Δ<0 | Empty |

At a double root, F(X)=−α(X−b/(2α))², so a nonconstant radial path is impossible.
The critical ray can still have angular motion before reduction.
At the singular orbit's axis crossing, (X,Y,Z)=(0,b,0); this is not the vertex
unless b=0. A small positive C is never classified as C=0.

With the page's positive Hamiltonian and Poisson brackets:

    dX/dz = −2Z/H
    dY/dz = 2αZ/H
    dZ/dz = (αX−Y)/H

The exact orbit uses θ=2√α z/H+φ, width=X₊−X₋:

    X = X₋ + width cos²(θ/2)
    Y = α(X₋ + width sin²(θ/2))
    Z = √α width sin(θ/2) cos(θ/2)
    T = πH/√α

Thus XY−Z²=αX₋X₊=C and αX+Y=b. No negative Y or F(X) is clipped to fabricate a
path. The double-root comparison uses 64 machine epsilons times
max(b²,4αC), and a second constraint check gates rendering. A negative Δ within
that floating-point bound is treated as the same numerically unresolved double
root; there is no 1e−5-wide physical "critical" band.

The critical button resets to n₀=1.5, α=.08, H=1.3 and computes C=b²/(4α)
without knob-step rounding. This avoids clamping an infeasible critical preset
at arbitrary previous knob settings. Ordinary knob interaction retains its
steps. Empty and stationary states disable playback while preserving the user's
play/pause preference for the next moving state.

## Run

Node 22 or newer; numerical tests require no packages:

    node --test tests/grin.test.cjs

Real Chromium tests:

    npm ci --prefix tests
    npx --prefix tests playwright install --with-deps chromium
    node tests/grin-browser.cjs

The browser script serves the production page on loopback HTTP. It exercises
1440×1000 desktop and 390×844 mobile/touch emulation, presets, the reported
C=.8/H=1.4 empty case, playback transitions, Tab/Enter, arrows, PageUp/PageDown,
Home/End, mouse/touch knobs, camera drag/reset, page errors and horizontal
overflow. It writes screenshots and browser-results.json to tests/artifacts/.
The GitHub Actions workflow runs both suites and uploads those artifacts.
Mobile emulation is not a physical phone, iOS Safari, or a screen-reader audit.

## Numerical comparison

Baseline: main commit 99b682e08647c71f9e1c0ebf18b33e08b5a7b924.
Values below use exact reported parameter inputs, rather than the rounded knob
defaults. Reproduce with:

    git show 99b682e08647c71f9e1c0ebf18b33e08b5a7b924:Web/poisson_variety_grin.html > /tmp/grin-baseline.html
    node tests/grin-compare.cjs /tmp/grin-baseline.html

All cases use α=.08, n₀=1.5, φ=0. Energy error is |n₀²−αX−Y−H²|.

| Case (C, H) | Points before → after | Max Casimir error before → after | Max energy error before → after |
| --- | --- | --- | --- |
| Guided (.030625, 1.4517) | 420 → 420 | 4.51e−17 → 2.08e−17 | 0 → 4.44e−16 |
| Critical (.98, 1.30) | 420 → 1 | .0157692 → 6.66e−16 | 0 → 2.22e−16 |
| Singular (0, 1.42) | 420 → 420 | 2.78e−17 → 8.33e−17 | 0 → 4.44e−16 |
| Empty (.8, 1.40) | 420 → 0 | .8 → N/A | .1116 → N/A |
| Singular, invalid energy (0, 1.58) | 420 → 0 | 0 → N/A | .36 → N/A |

N/A means no physical state exists; it is not a zero residual.

## Verification scope

The Node suite has 11 regression tests, including a 270-case parameter grid,
near-critical cases on both sides, all sampled energy curves, turning points,
the singular vertex, finite-difference checks of the reduced ODEs, an injected
invalid state to exercise the rendering guard, and animation/keyboard handlers.
Its DOM and canvas are doubles: this suite alone does not establish browser
rendering or visual correctness.

Use the CI job and browser-results.json to determine whether the real browser
suite ran successfully. A missing browser executable or a browser launch failure
must be reported as an unexecuted browser test, not as a pass.

Remaining model limits: the UI supports α>0, H>0 and C≥0 only. The plotted 420
samples cover one reduced radial period; straight canvas segments approximate
the continuous curve. The visualization does not reconstruct an unreduced
three-dimensional optical ray.
