# Analytical solution code

This folder contains the MATLAB code that evaluates the analytical
solution for the mean eccentricity vector (Cook's zonal framework
extended with odd-degree tesseral/sectorial forcing terms, as derived in
the paper), together with the auxiliary orbit-mechanics and special-
function routines it depends on.

## Requirements

- MATLAB. No toolboxes beyond base MATLAB are required (Parallel
  Computing Toolbox is only used for the `parfor` loop in
  `run_analytical_solution.m`; replace it with a plain `for` if not
  available — the loop is embarrassingly parallel either way).

## Layout

- `matlab/run_analytical_solution.m` — evaluates the analytical solution
  over a grid of initial (eccentricity, inclination) conditions at a fixed
  altitude and saves the resulting (a, e, i, RAAN, AOP) time histories.
  This mirrors the numerical grid campaign in
  `../propagation/matlab/propagate_grid_campaign.m` and produced the
  analytical side of the eccentricity-error maps in the paper.
- `matlab/AnalyticalSolution.m` — the analytical solution itself: given
  mean Keplerian elements at epoch, returns e, AOP, RAAN, i, a as a
  function of time.
- `matlab/CoefMoon.m` — loads and denormalizes the lunar gravity
  coefficients (zonal Jn and the tesseral Clm/Slm used by the odd-degree
  forcing terms) from `data/Coef_moon.mat`.
- `matlab/mean2osc_KL_moon.m`, `KepEqn.m`, `Lyd_correction.m` — first-order
  Kozai/Lyddane mean-to-osculating element conversion (and its Kepler
  equation solver).
- `matlab/compute_Omega_dot.m`, `Kaula_INC_dot.m` — secular RAAN
  precession rate from the even zonal harmonics.
- `matlab/RACook_const.m`, `coeff_odd.m`, `Kaula_INC.m` — Cook's zonal
  frozen-eccentricity constants and the odd-tesseral forcing coefficients
  entering the eccentricity-vector equations.
- `matlab/legendre_rec.m`, `legendre_m1.m`, `legendre_m2.m`,
  `DerivadasLegendre.m` — associated Legendre polynomials/derivatives used
  by the functions above.
- `data/Coef_moon.mat` — lunar gravity coefficients (`Clm`, `Slm`) used by
  `CoefMoon.m`.

## Usage

Run `matlab/run_analytical_solution.m` from MATLAB with `matlab/` on the
path (or as the working directory).

## Notes for reviewers

- The original development pipeline ran a MATLAB-Coder-generated MEX file
  (`AnalyticalSolution_mex`, compiled for Windows, `.mexw64`) for speed on
  the full 100x100 grid. That binary is a build artifact, not source code,
  and is platform-specific, so it is not included here — this release
  calls the plain `AnalyticalSolution.m` function directly, which gives
  identical results, just slower on the full grid. If speed matters,
  regenerate the MEX file yourself with
  `codegen AnalyticalSolution -args {kep, tau, theta0, deg}` (MATLAB
  Coder required).
- Only the code path actually used by `run_analytical_solution.m` is
  included. Other exploratory/plotting scripts and alternate solution
  variants present in the authors' working repository are not part of
  this release.
