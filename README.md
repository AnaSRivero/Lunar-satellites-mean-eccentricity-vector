# Lunar eccentricity-vector code

Code accompanying "A simple analytical solution for the evolution of the
mean eccentricity vector in extremely low lunar orbits", provided for
reproducibility of the numerical and analytical results.

- [`propagation/`](propagation/) — Orekit-based numerical propagator used
  as the reference (osculating) solution.
- [`analytical_solution/`](analytical_solution/) — MATLAB implementation
  of the analytical solution.

Each folder has its own README with requirements and usage instructions.

## Reproducing an eccentricity-error map

For one altitude (25, 50 or 75 km), run in this order, setting the same
`delta_a` and file names in each script:

1. `propagation/matlab/propagate_grid_campaign.m` — numerical reference
   over the grid of initial conditions (365 days or until impact).
2. `propagation/matlab/propagate_initial_revolution_campaign.m` — the same
   dynamics over one revolution centred at the initial epoch, used to
   estimate the mean initial elements.
3. `analytical_solution/matlab/run_analytical_solution.m` — analytical
   solution over the same grid.
4. `analytical_solution/matlab/eccentricity_error_map.m` — maximum
   eccentricity difference between the two solutions and its map.
