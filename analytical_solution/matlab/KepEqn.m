function nu = KepEqn(M, e, units)
%KEPEQN Solve Kepler equation and return true anomaly.

if nargin < 3
    units = 'rad';
end

if strcmpi(units, 'deg')
    M = deg2rad(M);
end

M = mod(M, 2*pi);
E = M;

for kk = 1:50
    f = E - e.*sin(E) - M;
    fp = 1 - e.*cos(E);
    dE = -f./fp;
    E = E + dE;
    if all(abs(dE) < 1e-13)
        break
    end
end

nu = 2.*atan2(sqrt(1+e).*sin(E./2), sqrt(1-e).*cos(E./2));
nu = mod(nu, 2*pi);

if strcmpi(units, 'deg')
    nu = rad2deg(nu);
end
end
