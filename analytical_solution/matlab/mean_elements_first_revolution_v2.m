function [kep_m, T0, mean_vars] = mean_elements_first_revolution_v2(kep, a_IC)
% MEAN_ELEMENTS_FIRST_REVOLUTION_V2
%
% Estimates mean orbital elements at the initial epoch by averaging the
% osculating Keplerian elements of a numerical propagation over one
% orbital revolution centred at that epoch.
%
% Instead of averaging e, omega, i and Omega separately, the following
% nonsingular variables are averaged (together with a):
%
%   xi  = e*cos(omega)
%   eta = e*sin(omega)
%
%   p   = tan(i/2)*sin(Omega)
%   q   = tan(i/2)*cos(Omega)
%
% The averaging interval is exactly one initial Keplerian orbital period
% T0, from -T0/2 to T0/2 with respect to the epoch.
%
%
% INPUT
% -----
% kep : 5xN matrix containing osculating Keplerian elements
%
%       kep(1,:) = a       [km]
%       kep(2,:) = e       [-]
%       kep(3,:) = i       [rad]
%       kep(4,:) = RAAN    [rad]
%       kep(5,:) = AOP     [rad]
%
%       sampled at the times
%           t = unique(sort([-T0/2:60:T0/2, -T0/2, 0, T0/2]))   [s]
%       i.e. every 60 s from -T0/2, plus the epoch and the end point
%       (see propagation/matlab/propagate_initial_revolution_campaign.m).
%
% a_IC : initial osculating semi-major axis [km], which defines
%        T0 = 2*pi*sqrt(a_IC^3/mu)
%
%
% OUTPUT
% ------
% kep_m : 5x1 vector containing mean Keplerian elements
%
%         kep_m(1) = a_m
%         kep_m(2) = e_m
%         kep_m(3) = i_m
%         kep_m(4) = RAAN_m
%         kep_m(5) = AOP_m
%
% T0 : initial Keplerian orbital period [s]
%
% mean_vars : structure containing the averaged nonsingular variables
%
%             mean_vars.xi
%             mean_vars.eta
%             mean_vars.p
%             mean_vars.q
%
%
% Convention used:
%
%   xi  = e*cos(omega)
%   eta = e*sin(omega)
%
%   p   = tan(i/2)*sin(Omega)
%   q   = tan(i/2)*cos(Omega)

    %% -------------------------------------------------------------
    % Constants
    % --------------------------------------------------------------

    muMoon = 4902.8001;       % Moon gravitational parameter [km^3/s^2]


    %% -------------------------------------------------------------
    % Input checks
    % --------------------------------------------------------------

    if size(kep,1) ~= 5
        error(['kep must be a 5xN matrix containing ' ...
               '[a; e; i; RAAN; AOP].']);
    end


    if size(kep,2) < 2
        error('The propagation must contain at least two samples.');
    end


    %% -------------------------------------------------------------
    % Osculating elements
    % --------------------------------------------------------------

    a    = kep(1,:);
    e    = kep(2,:);
    inc  = kep(3,:);
    RAAN = kep(4,:);
    AOP  = kep(5,:);


    T0 = 2*pi*sqrt(a_IC^3/muMoon);

    %% -------------------------------------------------------------
    % Convert to nonsingular variables
    % -------------------------------------------------------------

    % Eccentricity variables
    xi  = e .* cos(AOP);
    eta = e .* sin(AOP);

    % Inclination variables
    p = tan(inc/2) .* sin(RAAN);
    q = tan(inc/2) .* cos(RAAN);


    %% ------------------------------------------------------------
    % Time vector
    %% ------------------------------------------------------------

    halfT = T0/2; 
    tiempos = -halfT:60:halfT;

    t_rev = unique(sort([tiempos, -halfT, 0, halfT]));
    %% -------------------------------------------------------------
    % Temporal averages
    % -------------------------------------------------------------

    a_m = trapz(t_rev, a) / T0;

    xi_m  = trapz(t_rev, xi)  / T0;
    eta_m = trapz(t_rev, eta) / T0;

    p_m   = trapz(t_rev, p)   / T0;
    q_m   = trapz(t_rev, q)   / T0;


    %% -------------------------------------------------------------
    % Recover mean eccentricity and argument of periapsis
    % -------------------------------------------------------------

    e_m = hypot(xi_m, eta_m);

    AOP_m = atan2(eta_m, xi_m);

    AOP_m = mod(AOP_m, 2*pi);


    %% -------------------------------------------------------------
    % Recover mean inclination and RAAN
    % -------------------------------------------------------------

    tan_i2 = hypot(p_m, q_m);

    i_m = 2*atan(tan_i2);

    RAAN_m = atan2(p_m, q_m);

    RAAN_m = mod(RAAN_m, 2*pi);


    %% -------------------------------------------------------------
    % Output
    % -------------------------------------------------------------

    kep_m = [a_m;
             e_m;
             i_m;
             RAAN_m;
             AOP_m];


    %% Optional diagnostic output

    mean_vars.xi  = xi_m;
    mean_vars.eta = eta_m;

    mean_vars.p   = p_m;
    mean_vars.q   = q_m;

end