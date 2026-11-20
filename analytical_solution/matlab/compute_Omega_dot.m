function C_comp = compute_Omega_dot(a_m, e_m, i_m, param, deg)
% COMPUTE_OMEGA_DOT  Secular precession rate of the RAAN (dOmega/dt)
%
%   C_comp = compute_Omega_dot(a_m, e_m, i_m, param, deg)
%
%   Computes the secular variation of the ascending node due to the even
%   zonal harmonics of the lunar gravity field (up to degree 'deg'), to
%   first order in the harmonic coefficients, evaluated at the mean
%   elements. The coefficients are obtained internally from CoefMoon.m.
%
%   The evolution of the RAAN is obtained as:
%       Omega(tau) = Omega_0 + C_comp * tau
%   where tau = t / sqrt(Re^3/mu) is the dimensionless time.
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
%   a_m     [km]       Mean semi-major axis (scalar)
%   e_m     [-]        Mean eccentricity (scalar)
%   i_m     [rad]      Mean inclination (scalar)
%   param   struct     Parameters of the central body:
%                        param.Re  [km]        Equatorial radius
%   deg     [-]        Maximum degree of the even zonal harmonics included
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
%   C_comp  [rad per unit of dimensionless time]
%                      Secular rate dOmega/dtau evaluated at the mean elements
%
% -------------------------------------------------------------------------
% REFERENCES
% -------------------------------------------------------------------------
%   - Kaula, W.M. (1966). "Theory of Satellite Geodesy." Blaisdell.
%   - Kozai, Y. (1959). "The motion of a close earth satellite."
%     AJ, 64, 367.
%
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
%   * The eccentricity dependence is truncated at order e^2.
%   * The expression is singular at i = 0 (undefined node).
%   * The odd zonal harmonics do not contribute to the first-order secular
%     precession of the RAAN and are not included.

% =========================================================
% Auxiliary variables
% =========================================================
n = sqrt(1 ./ (a_m/param.Re).^3);   % unperturbed mean motion

% =========================================================
% First-order secular contribution of the even zonal
% harmonics J2, J4, ..., J_deg
% =========================================================
[~,~,~,C_mtx,~] = CoefMoon;
C0 = [0;C_mtx(:,1)];

sum_even = 0;

for i = 2:2:deg
    s = i/2;
    F_dot = Kaula_INC_dot(2*s, 0,s, i_m);
    sum_even = sum_even + F_dot*C0(i)*(1 + i*(i+1)/4*e_m^2);
end

C_comp = n/sin(i_m)/sqrt(1-e_m^2)*sum_even;

end