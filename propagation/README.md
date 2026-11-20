# Numerical propagation code

This folder contains the numerical (reference) propagation code used to
validate the analytical solution presented in the paper. It propagates a
satellite around the Moon under the full gravitational dynamics (lunar
gravity field up to a given degree/order, plus Earth/Sun third-body
attraction) using [Orekit](https://www.orekit.org/) through MATLAB's Java
interface.

## Requirements

- MATLAB (tested with R2025b).
- Orekit and Hipparchus jars (developed and verified against Orekit
  13.1.2/13.1.4 and Hipparchus 4.0.2/4.0.3; available from
  https://www.orekit.org/ and https://www.hipparchus.org/, or via Maven
  Central). Orekit is distributed under the Apache License 2.0 — it is
  **not** redistributed here, only referenced as a dependency.
- The `orekit-data` physical data set plus the lunar gravity field file —
  see `data/README.md`.

## Layout

- `matlab/setup_orekit.m` — adds the Orekit/Hipparchus jars to the MATLAB
  Java classpath and points Orekit at the local `orekit-data` folder. Call
  this once per session (and once per `parfor` worker).
- `matlab/getMoonPropagator.m` — builds a numerical propagator around the
  Moon (lunar gravity field of a given degree/order + optional Earth/Sun
  third-body attraction), using the full GRGM1200B field (all harmonics up
  to the requested degree/order). This is the propagator used for the
  reference propagations.
- `matlab/propagate_single_case.m` — propagates one initial condition and
  plots the osculating elements (a, e, i, RAAN, AOP) vs. time. This is the
  fully osculating reference used against the analytical (mean-element)
  solution.
- `matlab/propagate_grid_campaign.m` — batch-propagates a grid of initial
  (eccentricity, inclination) conditions at a fixed altitude and saves the
  full time histories; this produced the eccentricity-error maps used to
  characterize the domain of validity of the analytical solution. Edit
  `delta_a`/`Name` at the top to reproduce each altitude case (25/50/75 km).
- `matlab/propagate_initial_revolution_campaign.m` — for the same grid,
  force model and epoch, propagates one revolution centred at the initial
  epoch and saves the osculating elements every 60 s. These short
  histories are averaged by the analytical-solution code to estimate the
  mean initial elements, so this script must be run before
  `../analytical_solution/matlab/run_analytical_solution.m`.

## Usage

1. Edit the paths at the top of `propagate_single_case.m` /
   `propagate_grid_campaign.m` / `propagate_initial_revolution_campaign.m`
   (`orekitJarDir`, `orekitDataDir`) to point to your local Orekit jars
   and data folder.
2. Run the script from MATLAB with `matlab/` on the path.
