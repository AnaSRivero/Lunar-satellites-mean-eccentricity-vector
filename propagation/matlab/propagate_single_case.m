%% Single-orbit numerical propagation around the Moon
% Reproduces a single trajectory (a, e, i, RAAN, AOP vs. time) using the
% full gravitational dynamics: lunar gravity field up to (degree, order)
% plus Earth/Sun third-body attraction. This is the "reference"
% propagation used to validate the analytical solution in the paper.

clear; clc;

%% 1. Configuration (EDIT THESE PATHS)
orekitJarDir  = 'C:\Program Files\MATLAB\R2025b\java\jarext\';  % folder with orekit/hipparchus jars
orekitDataDir = 'C:\orekit-data-luna';                          % orekit-data + lunar_grail.gfc

setup_orekit(orekitJarDir, orekitDataDir);

import org.orekit.bodies.CelestialBodyFactory
import org.orekit.orbits.KeplerianOrbit
import org.orekit.orbits.PositionAngleType
import org.orekit.time.AbsoluteDate
import org.orekit.time.TimeScalesFactory
import org.orekit.utils.Constants

%% 2. Propagation window and perturbation model
initialDate = AbsoluteDate(2023, 1, 1, 12, 0, 0.0, TimeScalesFactory.getUTC());
duration    = 86400.0 * 200;          % [s]
finalDate   = initialDate.shiftedBy(duration);

degree    = 50;
order     = 50;
ThirdBody = 1;                        % 1: include Earth + Sun

%% 3. Initial orbit (Keplerian elements, Moon-centered inertial frame)
moon = CelestialBodyFactory.getMoon();
muMoon = moon.getGM();
moonInertialFrame = moon.getInertiallyOrientedFrame();

R_moon = Constants.MOON_EQUATORIAL_RADIUS;   % 1737400.0 m
a     = R_moon + 200e3;   % semi-major axis [m]
e     = 0;                % eccentricity
i     = deg2rad(27);      % inclination [rad]
omega = 0;                % argument of perigee [rad]
raan  = 0;                % RAAN [rad]
lv    = 0;                % true anomaly [rad]

initialOrbit = KeplerianOrbit(a, e, i, omega, raan, lv, ...
    PositionAngleType.TRUE, moonInertialFrame, initialDate, muMoon);

%% 4. Propagate
propagator = getMoonPropagator(initialOrbit, degree, order, ThirdBody);
generator  = propagator.getEphemerisGenerator();
propagator.propagate(initialDate, finalDate);
ephemeris  = generator.getGeneratedEphemeris();

%% 5. Sample the ephemeris
step_s = 600; % 10 min
t = 0:step_s:duration;
n = length(t);

res_a    = zeros(1, n);
res_e    = zeros(1, n);
res_i    = zeros(1, n);
res_AOP  = zeros(1, n);
res_RAAN = zeros(1, n);

for k = 1:n
    date_k  = initialDate.shiftedBy(t(k));
    state_k = ephemeris.propagate(date_k);
    orb_kep = org.orekit.orbits.KeplerianOrbit(state_k.getOrbit());

    res_a(k)    = orb_kep.getA() / 1000;   % km
    res_e(k)    = orb_kep.getE();
    res_i(k)    = rad2deg(orb_kep.getI());
    res_AOP(k)  = rad2deg(orb_kep.getPerigeeArgument());
    res_RAAN(k) = rad2deg(orb_kep.getRightAscensionOfAscendingNode());
end

%% 6. Plots
figure; plot(t/3600/24, res_e);    ylabel('Eccentricity [-]');       xlabel('Time [days]'); grid on
figure; plot(t/3600/24, res_i);    ylabel('Inclination [deg]');      xlabel('Time [days]'); grid on
figure; plot(t/3600/24, res_a);    ylabel('Semimajor axis [km]');    xlabel('Time [days]'); grid on
figure; plot(t/3600/24, res_RAAN); ylabel('RAAN [deg]');             xlabel('Time [days]'); grid on
figure; plot(t/3600/24, res_AOP);  ylabel('AOP [deg]');              xlabel('Time [days]'); grid on
