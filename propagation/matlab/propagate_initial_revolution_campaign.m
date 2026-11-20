%% One-revolution propagations for the mean initial elements
% For each initial condition (e0, i0) of the grid at a fixed altitude,
% propagates the reference dynamics over one Keplerian period centred at
% the initial epoch t0 (from t0 - T0/2 to t0 + T0/2) and saves the
% osculating (a, e, i, RAAN, AOP) histories sampled every 60 s.
%
% These short histories are the input of the analytical solution: they are
% averaged in analytical_solution/matlab/mean_elements_first_revolution_v2.m
% to estimate the mean initial elements (Section 5 of the paper). Run this
% script before analytical_solution/matlab/run_analytical_solution.m.
%
% Same force model, epoch and grid as propagate_grid_campaign.m: change
% delta_a / Name below to reproduce each altitude case (25, 50, 75 km).

clear; clc;

%% 1. Configuration (EDIT THESE PATHS)
orekitJarDir  = 'C:\Program Files\MATLAB\R2025b\java\jarext\';
orekitDataDir = 'C:\orekit-data-luna';
outputDir     = './Propagations';
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

%% 2. Campaign parameters
Name    = 'R75_IC';   % run label, used in the saved .mat file names
delta_a = 75;         % altitude above the lunar surface [km]
degree  = 50;
order   = 50;
ThirdBody = 1;        % 1: include Earth + Sun

R_moon = 1737400.0;   % [m]
ef     = 0.7 * (delta_a * 1e3 / (R_moon + delta_a * 1e3));  % 70% of impact eccentricity
grid   = 100;

a     = R_moon + delta_a * 1e3;         % semi-major axis [m]
e_v   = linspace(0, ef, grid);          % eccentricity grid
i_v   = deg2rad(linspace(0, 90, grid)); % inclination grid [rad]
omega = 0;
raan  = 0;
lv    = 0;

% One initial Keplerian period centred at the epoch
muMoon_T0 = 4902800118457.55;           % [m^3/s^2]
T0    = 2 * pi * sqrt(a^3 / muMoon_T0); % [s]
halfT = T0 / 2;

% Sample times relative to the epoch. This is the time vector assumed by
% mean_elements_first_revolution_v2.m: do not change one without the other.
step_s = 60;
t = -halfT:step_s:halfT;
t = unique(sort([t, -halfT, 0, halfT]));
n_steps = length(t);

e_prop    = zeros(length(e_v), n_steps, length(i_v));
AOP_prop  = zeros(length(e_v), n_steps, length(i_v));
a_prop    = zeros(length(e_v), n_steps, length(i_v));
i_prop    = zeros(length(e_v), n_steps, length(i_v));
RAAN_prop = zeros(length(e_v), n_steps, length(i_v));

%% 3. Run the grid (parallel over inclination for each eccentricity row)
% parpool('local', 8); % uncomment to control the number of workers

for ii = 1:length(e_v)
    e = e_v(ii);
    parfor jj = 1:length(i_v)
        setup_orekit(orekitJarDir, orekitDataDir); % each worker needs its own classpath

        import org.orekit.bodies.CelestialBodyFactory
        import org.orekit.orbits.KeplerianOrbit
        import org.orekit.orbits.PositionAngleType
        import org.orekit.time.AbsoluteDate
        import org.orekit.time.TimeScalesFactory

        initialDate = AbsoluteDate(2000, 1, 1, 12, 0, 0.0, TimeScalesFactory.getUTC());
        startDate   = initialDate.shiftedBy(-halfT);
        finalDate   = initialDate.shiftedBy(halfT);

        moon = CelestialBodyFactory.getMoon();
        muMoon = moon.getGM();
        moonInertialFrame = moon.getInertiallyOrientedFrame();

        i = i_v(jj);
        initialOrbit = KeplerianOrbit(a, e, i, omega, raan, lv, ...
            PositionAngleType.TRUE, moonInertialFrame, initialDate, muMoon);

        % The initial state is given at t0: the propagator first goes back
        % to t0 - T0/2 and then stores the ephemeris up to t0 + T0/2.
        propagator = getMoonPropagator(initialOrbit, degree, order, ThirdBody);
        generator  = propagator.getEphemerisGenerator();
        propagator.propagate(startDate, finalDate);
        ephemeris  = generator.getGeneratedEphemeris();

        res_e    = nan(1, n_steps);
        res_AOP  = nan(1, n_steps);
        res_a    = nan(1, n_steps);
        res_i    = nan(1, n_steps);
        res_RAAN = nan(1, n_steps);

        firstDate = ephemeris.getMinDate();
        lastDate  = ephemeris.getMaxDate();

        for k = 1:n_steps
            date_k = initialDate.shiftedBy(t(k));
            if date_k.compareTo(firstDate) >= 0 && date_k.compareTo(lastDate) <= 0
                state_k = ephemeris.propagate(date_k);
                orb_kep = org.orekit.orbits.KeplerianOrbit(state_k.getOrbit());

                res_AOP(k)  = rad2deg(orb_kep.getPerigeeArgument());
                res_e(k)    = orb_kep.getE();
                res_a(k)    = orb_kep.getA();
                res_i(k)    = rad2deg(orb_kep.getI());
                res_RAAN(k) = rad2deg(orb_kep.getRightAscensionOfAscendingNode());
            end
        end

        e_prop(ii, :, jj)    = res_e;
        AOP_prop(ii, :, jj)  = res_AOP;
        a_prop(ii, :, jj)    = res_a;
        i_prop(ii, :, jj)    = res_i;
        RAAN_prop(ii, :, jj) = res_RAAN;
    end
end

%% 4. Save results (a in m, angles in deg, as in propagate_grid_campaign.m)
save(fullfile(outputDir, ['Results_' Name]), ...
    'e_prop', 'AOP_prop', 'a_prop', 'i_prop', 'RAAN_prop', '-v7.3');
save(fullfile(outputDir, ['Configuration_' Name]), ...
    'a', 'e_v', 'i_v', 'T0', 't');
