%% Analytical solution: eccentricity/inclination grid
% Evaluates the analytical solution for the mean eccentricity vector (and
% the associated a, i, RAAN histories) over a grid of initial
% (eccentricity, inclination) conditions at a fixed altitude, matching the
% numerical grid campaign in ../../propagation/matlab/propagate_grid_campaign.m.
% This produced the analytical side of the eccentricity-error maps used to
% assess the domain of validity of the model (Section 6 of the paper).
%
% Calls AnalyticalSolution.m directly. The original development pipeline
% used a MATLAB-Coder-generated MEX file (AnalyticalSolution_mex, compiled
% for Windows) for speed on the full 100x100 grid; that binary is not
% portable and is not needed to reproduce the results, only to run them
% faster. Regenerate it yourself if needed with:
%   codegen AnalyticalSolution -args {kep, tau, theta0, deg}

clear; clc;

%% 1. Configuration (EDIT THIS PATH)
outputDir = './Results';
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

%% 2. Campaign parameters (mirrors the R75 case in the paper)
R_moon  = 1737.4;   % Moon's radius [km]
delta_a = 75;       % altitude above the lunar surface [km]
ef      = 0.7 * (delta_a / (R_moon + delta_a));  % 70% of impact eccentricity
grid    = 100;

a    = R_moon + delta_a;                 % semi-major axis [km]
e_v  = linspace(0, ef, grid);            % eccentricity grid
e_v(1) = 1e-10;                          % avoid the e=0 singularity in AOP
i_v  = deg2rad(linspace(0, 90, grid));   % inclination grid [rad]
i_v(1) = 1e-10;                          % avoid the i=0 singularity
AOP  = 0;   % argument of perigee [rad]
RAAN = 0;   % RAAN [rad]
MA   = 0;   % mean anomaly [rad]

deg = 50;   % max degree of the odd-tesseral/zonal series in the solution

param.Re = R_moon;                    % Moon's radius [km]
param.mu = 4.9028000661637961E+03;    % Moon's GM [km^3/s^2] (SPICE)
theta0   = 0.719163820108406;         % Moon's rotation angle at t0 (2000-01-01 12:00:00)
days     = 365;

duration = days * 24 * 3600;          % [s]
tev = 0:600:duration;                 % sample times [s], 10 min step
tau = (tev - tev(1)) ./ sqrt(param.Re^3 / param.mu); % dimensionless time

%% 3. Evaluate the analytical solution over the grid
e_analy_v    = zeros(length(e_v), length(tev), length(i_v));
AOP_analy_v  = zeros(length(e_v), length(tev), length(i_v));
RAAN_analy_v = zeros(length(e_v), length(tev), length(i_v));
i_analy_v    = zeros(length(e_v), length(tev), length(i_v));
a_analy_v    = zeros(length(e_v), length(tev), length(i_v));

for ii = 1:length(e_v)
    e = e_v(ii);
    parfor jj = 1:length(i_v)
        INC = i_v(jj);
        kep = [a'; e'; INC'; RAAN'; AOP'; MA'];

        [e_analy, AOP_analy, RAAN_analy, i_analy, a_analy] = AnalyticalSolution(kep, tau, theta0, deg, R_moon+delta_a);

        e_analy_v(ii, :, jj)    = e_analy;
        AOP_analy_v(ii, :, jj)  = AOP_analy;
        RAAN_analy_v(ii, :, jj) = RAAN_analy;
        i_analy_v(ii, :, jj)    = i_analy;
        a_analy_v(ii, :, jj)    = a_analy;
    end
end

%% 4. Save results
save(fullfile(outputDir, 'Results_R75'), ...
    'e_analy_v', 'AOP_analy_v', 'RAAN_analy_v', 'i_analy_v', 'a_analy_v', 'e_v', 'i_v', '-v7.3');
