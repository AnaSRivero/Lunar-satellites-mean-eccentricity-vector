function Pn2 = legendre_m2(n, x)
%#codegen
% Calcula P_n^2(x)
% Requiere n >= 2
Pn2 = 0.0; % Inicialización obligatoria

if n < 2
    Pn2 = 0.0;
    return
end

% P2^2
P2 = 3*(1 - x*x);

if n == 2
    Pn2 = P2;
    return
end

% P3^2
P3 = 5*x*P2;

if n == 3
    Pn2 = P3;
    return
end

Pnm2 = P2;
Pnm1 = P3;

for k = 4:n
    Pn2 = ((2*k-1)*x*Pnm1 - (k+1)*Pnm2)/(k-2);
    Pnm2 = Pnm1;
    Pnm1 = Pn2;
end
