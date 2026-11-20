%% Analytical solution: eccentricity/inclination grid
% Evaluates the analytical solution for the mean eccentricity vector (and
% the associated a, i, RAAN histories) over the grid of initial
% (eccentricity, inclination) conditions at a fixed altitude used in the
% numerical grid campaign (../../propagation/matlab/propagate_grid_campaign.m).
% This produces the analytical side of the eccentricity-error maps used to
% assess the domain of applicability of the model (Section 5 of the paper).
%
% The mean initial elements are estimated numerically: for each grid point,
% the osculating elements of the reference propagation over one revolution
% centred at the initial epoch are averaged inside AnalyticalSolution.m
% (see mean_elements_first_revolution_v2.m). Those one-revolution histories
% must be generated first with
%   ../../propagation/matlab/propagate_initial_revolution_campaign.m
%
% Calls AnalyticalSolution.m directly. The original development pipeline
% used a MATLAB-Coder-generated MEX file for speed on the full 100x100
% grid; that binary is not portable and is not needed to reproduce the
% results, only to run them faster.

clear; clc;

%% 1. Configuration (EDIT THESE PATHS)
icFile    = '../../propagation/matlab/Propagations/Results_R75_IC.mat'; % one-revolution histories
outputDir = './Results';
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

%% 2. Campaign parameters (must match the run that produced icFile)
R_moon  = 1737.4;   % Moon's radius [km]
delta_a = 75;       % altitude above the lunar surface [km]
ef      = 0.7 * (delta_a / (R_moon + delta_a));  % 70% of impact eccentricity
grid    = 100;

a_IC = R_moon + delta_a;                 % initial osculating semi-major axis [km]
e_v  = linspace(0, ef, grid);            % initial osculating eccentricity grid
i_v  = deg2rad(linspace(0, 90, grid));   % initial osculating inclination grid [rad]

deg = 50;   % max degree of the lunar gravity field in the solution

param.Re = R_moon;                    % Moon's radius [km]
param.mu = 4.9028000661637961E+03;    % Moon's GM [km^3/s^2] (SPICE)
theta0   = 0.719163820108406;         % Moon's rotation angle at t0 (2000-01-01 12:00:00)
days     = 365;

duration = days * 24 * 3600;          % [s]
tev = 0:600:duration;                 % sample times [s], 10 min step
tau = (tev - tev(1)) ./ sqrt(param.Re^3 / param.mu); % dimensionless time

%% 3. One-revolution osculating histories (a [m], angles [deg])
IC = load(icFile, 'a_prop', 'e_prop', 'i_prop', 'RAAN_prop', 'AOP_prop');
if size(IC.e_prop, 1) ~= length(e_v) || size(IC.e_prop, 3) ~= length(i_v)
    error('run_analytical_solution:gridMismatch', ...
        'The grid in %s does not match the campaign parameters above.', icFile);
end

%% 4. Evaluate the analytical solution over the grid
e_analy_v    = zeros(length(e_v), length(tev), length(i_v));
AOP_analy_v  = zeros(length(e_v), length(tev), length(i_v));
RAAN_analy_v = zeros(length(e_v), length(tev), length(i_v));
i_analy_v    = zeros(length(e_v), length(tev), length(i_v));
a_analy_v    = zeros(length(e_v), length(tev), length(i_v));

for ii = 1:length(e_v)
    % histories of this eccentricity row: (samples) x (inclinations)
    a_row    = squeeze(IC.a_prop(ii, :, :)) ./ 1000;      % [km]
    e_row    = squeeze(IC.e_prop(ii, :, :));
    i_row    = deg2rad(squeeze(IC.i_prop(ii, :, :)));
    RAAN_row = deg2rad(squeeze(IC.RAAN_prop(ii, :, :)));
    AOP_row  = deg2rad(squeeze(IC.AOP_prop(ii, :, :)));

    parfor jj = 1:length(i_v)
        kep = [a_row(:, jj)'; e_row(:, jj)'; i_row(:, jj)'; RAAN_row(:, jj)'; AOP_row(:, jj)'];

        [e_analy, AOP_analy, RAAN_analy, i_analy, a_analy] = AnalyticalSolution(kep, tau, theta0, deg, a_IC);

        e_analy_v(ii, :, jj)    = e_analy;
        AOP_analy_v(ii, :, jj)  = AOP_analy;
        RAAN_analy_v(ii, :, jj) = RAAN_analy;
        i_analy_v(ii, :, jj)    = i_analy;
        a_analy_v(ii, :, jj)    = a_analy;
    end
end

%% 5. Save results
save(fullfile(outputDir, 'Results_R75'), ...
    'e_analy_v', 'AOP_analy_v', 'RAAN_analy_v', 'i_analy_v', 'a_analy_v', 'e_v', 'i_v', '-v7.3');
