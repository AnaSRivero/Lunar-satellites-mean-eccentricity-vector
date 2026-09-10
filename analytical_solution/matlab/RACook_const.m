function [k_xi, k_eta, C] = RACook_const (SMA_KM,INC_DEG,N_H)
%--------------------------------------------------------------------------
%  Computes Cook's zonal eccentricity-vector constants k_xi, k_eta (the
%  homogeneous-solution frequencies, whose sign/product also determines
%  the frozen-eccentricity behaviour) and the odd-zonal coupling term C,
%  from the mean semi-major axis and inclination, considering N_H zonal
%  harmonics of the lunar gravity field.
%--------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   SMA_KM    (1,1)   mean semimajor axis [km]
%   INC_DEG   (1,1)   mean inclination [deg]
%   N_H       (1,1)   number of zonal harmonics to include (max 23,
%                      limited by the coefficients in CoefMoon.m)
%
%   -------
%   Outputs
%   -------
%   k_xi, k_eta (1,1)  eccentricity-vector frequencies [1/s]
%   C           (1,1)  odd-zonal coupling term [1/s]
%
%--------------------------------------------------------------------------


% Moon radius (km)
RE= 1737.4;               % Moon's radius [km]

[J,~,~] = CoefMoon;
a = SMA_KM;
n = sqrt(1./(a./RE).^3);
INC = INC_DEG/180*pi;


% Computing Cook's frecuency of xi and eta
k_xi = 0;
k_eta = 0;
for ii = 2:2:N_H
    Pn_cos = legendre_rec(ii, cos(INC));
    Pn_0 = legendre_rec(ii,0);
    Pn2_cos = legendre_m2(ii, cos(INC));
    Pn2_0 = legendre_m2(ii,0);
    Pnd = DerivadasLegendre(INC, ii);
    k_xi = k_xi + J(ii)*n*(RE/a)^ii*(((ii-1)*(ii-2)*factorial(ii-2))/(2*factorial(ii+2))*Pn2_cos*Pn2_0 + ii*(ii+1)/2*Pn_cos*Pn_0 - Pnd*Pn_0*(1/tan(INC)));
    k_eta = k_eta + J(ii)*n*(RE/a)^ii*(((ii-1)*(ii-2)*factorial(ii-2))/(2*factorial(ii+2))*Pn2_cos*Pn2_0 - ii*(ii+1)/2*Pn_cos*Pn_0 + Pnd*Pn_0*(1/tan(INC)));
end

C = 0;
for ii = 3:2:N_H
    Pn1_cos = legendre_m1(ii,cos(INC));
    Pn1_0 = legendre_m1(ii,0);

    C = C + J(ii).*n.*(RE./a).^ii.*(ii-1)./(ii.*(ii+1)).*Pn1_cos.*Pn1_0;
end

end
