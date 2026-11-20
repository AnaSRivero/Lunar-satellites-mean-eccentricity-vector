%% Eccentricity-error map
% Compares the analytical solution with the numerical reference over the
% grid of initial conditions and plots the validation metric of the paper
% (Section 5):
%
%   Delta e_max = max over t of |e_analytical(t) - e_numerical(t)|
%
% for each initial (eccentricity, inclination), with t from the initial
% epoch to 365 days or to the impact time of the numerical trajectory,
% whichever occurs first. The map shows log10(Delta e_max).
%
% Inputs (run these first, with the same altitude):
%   ../../propagation/matlab/propagate_grid_campaign.m   -> numerical reference
%   run_analytical_solution.m                            -> analytical solution

clear; clc;

%% 1. Configuration (EDIT THESE PATHS)
numericalFile  = '../../propagation/matlab/Propagations/Results_R75.mat';
analyticalFile = './Results/Results_R75.mat';

R_moon  = 1737.4;   % Moon's radius [km]
delta_a = 75;       % altitude above the lunar surface [km]
ef      = 0.7 * (delta_a / (R_moon + delta_a));
grid    = 100;

e_v = linspace(0, ef, grid);    % initial osculating eccentricity
i_v = linspace(0, 90, grid);    % initial osculating inclination [deg]

%% 2. Maximum eccentricity difference over the comparison interval
N = load(numericalFile,  'e_prop');       % (eccentricity, time, inclination)
A = load(analyticalFile, 'e_analy_v');    % (eccentricity, time, inclination)

n_t = min(size(N.e_prop, 2), size(A.e_analy_v, 2));
e_num = N.e_prop(:, 1:n_t, :);
e_num(e_num <= 0) = NaN;                  % samples after impact are not compared

delta_e_max = squeeze(max(abs(e_num - A.e_analy_v(:, 1:n_t, :)), [], 2, 'omitnan'));

%% 3. Map
figure;
imagesc(i_v, e_v, log10(delta_e_max));
axis xy; axis square;
clim([-3 -1]);
h = colorbar;
ylabel(h, 'max eccentricity error (log)');
xlabel('Inclination (deg)');
ylabel('Eccentricity (-)');
xlim([0 90]); ylim([0 max(e_v)]);
