%% Eccentricity/inclination error-map campaign
% Propagates a grid of initial conditions (e0, i0) at a fixed altitude and
% saves the full (a, e, i, RAAN, AOP) time histories. This is the batch
% script behind the eccentricity-error maps used to characterize the
% domain of validity of the analytical solution (Section 6 of the paper).
%
% Generalizes the original 'Campain_propagation.m' / '_server.m' scripts:
% just change delta_a / Name below to reproduce each altitude case
% (e.g. 25, 50, 75 km).

clear; clc;

%% 1. Configuration (EDIT THESE PATHS)
orekitJarDir  = 'C:\Program Files\MATLAB\R2025b\java\jarext\';
orekitDataDir = 'C:\orekit-data-luna';
outputDir     = './Propagations';
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

%% 2. Campaign parameters
Name    = 'R50';      % run label, used in the saved .mat file names
delta_a = 50;         % altitude above the lunar surface [km]
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

days     = 365;
duration = 86400.0 * days;

step_s   = 600; % 10 min
n_steps  = floor(duration / step_s) + 1;

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
        finalDate   = initialDate.shiftedBy(duration);

        moon = CelestialBodyFactory.getMoon();
        muMoon = moon.getGM();
        moonInertialFrame = moon.getInertiallyOrientedFrame();
        date0 = initialDate;

        i = i_v(jj);
        initialOrbit = KeplerianOrbit(a, e, i, omega, raan, lv, ...
            PositionAngleType.TRUE, moonInertialFrame, date0, muMoon);

        propagator = getMoonPropagator(initialOrbit, degree, order, ThirdBody);
        generator  = propagator.getEphemerisGenerator();
        try
            propagator.propagate(initialDate, finalDate);
        catch ME
            if contains(ME.message, 'minimal step size')
                warning('Impact detected: propagation stopped early.');
            else
                rethrow(ME);
            end
        end
        ephemeris = generator.getGeneratedEphemeris();

        t = 0:step_s:duration;
        n = length(t);

        res_e    = nan(1, n);
        res_AOP  = nan(1, n);
        res_a    = nan(1, n);
        res_i    = nan(1, n);
        res_RAAN = nan(1, n);

        firstDate = ephemeris.getMinDate();
        lastDate  = ephemeris.getMaxDate();

        for k = 1:n
            date_k = date0.shiftedBy(t(k));
            if date_k.compareTo(firstDate) >= 0 && date_k.compareTo(lastDate) <= 0
                state_k = ephemeris.propagate(date_k);
                orb_kep = org.orekit.orbits.KeplerianOrbit(state_k.getOrbit());

                res_AOP(k)  = rad2deg(orb_kep.getPerigeeArgument());
                res_e(k)    = orb_kep.getE();
                res_a(k)    = orb_kep.getA();
                res_i(k)    = rad2deg(orb_kep.getI());
                res_RAAN(k) = rad2deg(orb_kep.getRightAscensionOfAscendingNode());

                pos = state_k.getPVCoordinates().getPosition().getNorm();
                if pos - R_moon <= 0
                    break; % impact: remaining samples stay NaN
                end
            end
        end

        e_prop(ii, :, jj)    = res_e;
        AOP_prop(ii, :, jj)  = res_AOP;
        a_prop(ii, :, jj)    = res_a;
        i_prop(ii, :, jj)    = res_i;
        RAAN_prop(ii, :, jj) = res_RAAN;
    end
end

%% 4. Save results
save(fullfile(outputDir, ['Results_' Name]), ...
    'e_prop', 'AOP_prop', 'a_prop', 'i_prop', 'RAAN_prop', '-v7.3');
save(fullfile(outputDir, ['Configuration_' Name]), ...
    'a', 'e_v', 'i_v', 'duration');
