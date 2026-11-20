%% Reduced-model propagations: effect of omitted even-degree harmonics
% Propagates the two cases of Section 5 of the paper in which the
% analytical solution departs from the numerical reference, using reduced
% numerical models of the lunar gravity field:
%
%   'reduced'   : all zonal harmonics + odd-degree non-zonal harmonics
%                 (the harmonics retained in the analytical solution)
%   'augmented' : the same model plus selected even-degree non-zonal
%                 harmonics (C_lm and S_lm)
%
% Case 1: 50 km altitude, e0 = 0.005736, I0 = 49.09 deg; augmented with
%         C22, S22, C43, S43.
% Case 2: 75 km altitude, e0 = 0.008485, I0 = 73.64 deg; augmented with
%         C22, S22.
%
% Both are propagated for 365 days up to degree and order 50, without
% third-body attraction, and the osculating (a, e, i, RAAN, AOP) are saved
% every 600 s. The full reference trajectories of the same initial
% conditions are obtained with propagate_single_case.m or from the grid
% campaign.

clear; clc;

%% 1. Configuration (EDIT THESE PATHS)
orekitJarDir  = 'C:\Program Files\MATLAB\R2025b\java\jarext\';
orekitDataDir = 'C:\orekit-data-luna';
outputDir     = './Propagations';
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

setup_orekit(orekitJarDir, orekitDataDir);

import org.orekit.bodies.CelestialBodyFactory
import org.orekit.orbits.KeplerianOrbit
import org.orekit.orbits.PositionAngleType
import org.orekit.time.AbsoluteDate
import org.orekit.time.TimeScalesFactory

%% 2. Cases and models
R_moon = 1737400.0;   % [m]
degree = 50;
order  = 50;
days   = 365;
grid   = 100;         % grid of Section 5, used to define case 1

i_v = linspace(0, 90, grid);                       % [deg]

h = 50;  e_v = linspace(0, 0.7 * (h * 1e3 / (R_moon + h * 1e3)), grid);
cases(1).name = 'case1_50km';
cases(1).h = h;  cases(1).e = e_v(30);  cases(1).i = i_v(55);
cases(1).extraTerms = [2 2; 4 3];

cases(2).name = 'case2_75km';
cases(2).h = 75;  cases(2).e = 0.008485;  cases(2).i = 73.64;
cases(2).extraTerms = [2 2];

duration = 86400 * days;      % [s]
t = 0:600:duration;           % sample times [s]
n = numel(t);

initialDate = AbsoluteDate(2000, 1, 1, 12, 0, 0.0, TimeScalesFactory.getUTC());
moon   = CelestialBodyFactory.getMoon();
muMoon = moon.getGM();
moonInertialFrame = moon.getInertiallyOrientedFrame();

%% 3. Propagate
res = struct();
for c = 1:numel(cases)
    models = {'reduced', zeros(0, 2); 'augmented', cases(c).extraTerms};
    for m = 1:size(models, 1)
        a = R_moon + cases(c).h * 1e3;
        initialOrbit = KeplerianOrbit(a, cases(c).e, deg2rad(cases(c).i), 0, 0, 0, ...
            PositionAngleType.TRUE, moonInertialFrame, initialDate, muMoon);

        propagator = getMoonPropagatorReduced(initialOrbit, degree, order, models{m, 2}, orekitDataDir);
        generator  = propagator.getEphemerisGenerator();
        try
            propagator.propagate(initialDate, initialDate.shiftedBy(duration));
        catch ME
            % the integration stops at impact with a minimal-step-size error
            if ~contains(ME.message, 'minimal step size'), rethrow(ME); end
        end
        ephemeris = generator.getGeneratedEphemeris();

        out = struct('t', t, 'a', nan(1, n), 'e', nan(1, n), 'i', nan(1, n), ...
            'RAAN', nan(1, n), 'AOP', nan(1, n));
        for k = 1:n
            date_k = initialDate.shiftedBy(t(k));
            if date_k.compareTo(ephemeris.getMinDate()) >= 0 && date_k.compareTo(ephemeris.getMaxDate()) <= 0
                state_k = ephemeris.propagate(date_k);
                orb_kep = KeplerianOrbit(state_k.getOrbit());
                out.a(k)    = orb_kep.getA();                                        % [m]
                out.e(k)    = orb_kep.getE();
                out.i(k)    = rad2deg(orb_kep.getI());                               % [deg]
                out.RAAN(k) = rad2deg(orb_kep.getRightAscensionOfAscendingNode());   % [deg]
                out.AOP(k)  = rad2deg(orb_kep.getPerigeeArgument());                 % [deg]
            end
        end
        res.(cases(c).name).(models{m, 1}) = out;
        fprintf('%s / %s done\n', cases(c).name, models{m, 1});
        save(fullfile(outputDir, 'Results_reduced_models'), 'res', 'cases');
    end
end
