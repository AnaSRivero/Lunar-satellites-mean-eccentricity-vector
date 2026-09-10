function [QC_xi, QS_xi, QC_eta, QS_eta] = coeff_odd(a_m, i_m, n, beta_0, C_mtx, S_mtx, i, j)
%COEFF_ODD Forcing coefficients of the odd-degree tesseral/sectorial term
%(degree i, order j) entering the xi/eta eccentricity-vector equations.
%
%   [QC_xi, QS_xi, QC_eta, QS_eta] = coeff_odd(a_m, i_m, n, beta_0, ...
%       C_mtx, S_mtx, i, j)
%
%   a_m, i_m : mean semi-major axis [km] and inclination [rad]
%   n        : mean motion [1/s]
%   beta_0   : RAAN_m - theta0, longitude offset at epoch [rad]
%   C_mtx, S_mtx : denormalized tesseral coefficients (from CoefMoon.m)
%   i, j     : harmonic degree and order

param.Re = 1738;               % Moon's radius [km]

l = i-1;

if mod((i-j),2) ==  0
    C = C_mtx(l,j+1);
    S = S_mtx(l,j+1);
else
    C = -S_mtx(l,j+1);
    S = C_mtx(l,j+1);
end

F1 = Kaula_INC(i, j, (i-1)/2, i_m);
F3 = Kaula_INC(i, j, (i+1)/2, i_m);
gamma_n = (F1-F3)/2;
gamma_p = (F1+F3)/2;

m = j;

QC_xi = n.*(param.Re./a_m).^i.*(i - 1).*gamma_n.*(C*sin(m*beta_0)-S*cos(m*beta_0));
QS_xi = -n.*(param.Re./a_m).^i.*(i - 1).*gamma_n.*(C*cos(m*beta_0)+S*sin(m*beta_0));
QC_eta = n.*(param.Re./a_m).^i.*(i - 1).*gamma_p.*(C*cos(m*beta_0)+S*sin(m*beta_0));
QS_eta = n.*(param.Re./a_m).^i.*(i - 1).*gamma_p.*(C*sin(m*beta_0)-S*cos(m*beta_0));

end