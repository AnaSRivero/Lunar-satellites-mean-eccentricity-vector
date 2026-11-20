# Analytical solution code

This folder contains the MATLAB code that evaluates the analytical
solution for the mean eccentricity vector (Cook's zonal framework
extended with odd-degree tesseral/sectorial forcing terms, as derived in
the paper), together with the auxiliary special-function routines it
depends on.

## Requirements

- MATLAB. No toolboxes beyond base MATLAB are required (Parallel
  Computing Toolbox is only used for the `parfor` loop in
  `run_analytical_solution.m`; replace it with a plain `for` if not
  available — the loop is embarrassingly parallel either way).
- The one-revolution reference propagations produced by
  `../propagation/matlab/propagate_initial_revolution_campaign.m`, which
  require Orekit (see `../propagation/README.md`).

## Layout

- `matlab/run_analytical_solution.m` — evaluates the analytical solution
  over a grid of initial (eccentricity, inclination) conditions at a fixed
  altitude and saves the resulting (a, e, i, RAAN, AOP) time histories.
  This mirrors the numerical grid campaign in
  `../propagation/matlab/propagate_grid_campaign.m` and produces the
  analytical side of the eccentricity-error maps in the paper.
- `matlab/eccentricity_error_map.m` — computes the validation metric of
  the paper, the maximum absolute difference between the analytical and
  numerical eccentricities over the comparison interval, for every grid
  point, and plots its map.
- `matlab/AnalyticalSolution.m` — the analytical solution itself: given
  the osculating elements of the reference propagation over one
  revolution centred at the initial epoch, it estimates the mean initial
  elements and returns e, AOP, RAAN, i, a as a function of time.
- `matlab/mean_elements_first_revolution_v2.m` — estimates the mean
  initial elements by averaging the semi-major axis and the non-singular
  variables e cos(AOP), e sin(AOP), tan(i/2) sin(RAAN), tan(i/2) cos(RAAN)
  over that revolution.
- `matlab/CoefMoon.m` — loads and denormalizes the lunar gravity
  coefficients (zonal Jn and the tesseral Clm/Slm used by the odd-degree
  forcing terms) from `data/Coef_moon.mat`.
- `matlab/compute_Omega_dot.m`, `Kaula_INC_dot.m` — secular RAAN
  precession rate from the even zonal harmonics.
- `matlab/RACook_const.m`, `coeff_odd.m`, `Kaula_INC.m` — Cook's zonal
  frozen-eccentricity constants and the odd-tesseral forcing coefficients
  entering the eccentricity-vector equations. `Kaula_INC.m` evaluates the
  inclination functions through numerically stable closed forms;
  `Kaula_INC_series.m` is the direct series, kept as a fallback, which
  loses precision in double arithmetic for degrees above about 20.
- `matlab/legendre_rec.m`, `legendre_m1.m`, `legendre_m2.m`,
  `DerivadasLegendre.m` — associated Legendre polynomials/derivatives used
  by the functions above.
- `data/Coef_moon.mat` — lunar gravity coefficients (`Clm`, `Slm`) used by
  `CoefMoon.m`.

## Usage

1. Generate the one-revolution histories with
   `../propagation/matlab/propagate_initial_revolution_campaign.m`.
2. Set `icFile` and the campaign parameters at the top of
   `matlab/run_analytical_solution.m` and run it from MATLAB with
   `matlab/` as the working directory.
3. To obtain the error map, run `matlab/eccentricity_error_map.m` after
   the numerical grid campaign of the same altitude.

## Notes for reviewers

- The mean initial elements are numerical estimates obtained from the
  reference propagation, not an exact transformation to the mean elements
  of the analytical theory.
- The original development pipeline ran a MATLAB-Coder-generated MEX file
  (compiled for Windows, `.mexw64`) for speed on the full 100x100 grid.
  That binary is a build artifact, not source code, and is
  platform-specific, so it is not included here — this release calls the
  plain `AnalyticalSolution.m` function directly, which gives identical
  results, just slower on the full grid.
- Only the code path actually used by `run_analytical_solution.m` is
  included. Other exploratory/plotting scripts and alternate solution
  variants present in the authors' working repository are not part of
  this release.
