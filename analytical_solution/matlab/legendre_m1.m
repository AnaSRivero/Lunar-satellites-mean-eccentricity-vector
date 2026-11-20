function Pn1 = legendre_m1(n, x)
%#codegen
% Computes the associated Legendre function P_n^1(x)
Pn1 = 0.0; % mandatory initialization (code generation)

if n < 1
    Pn1 = 0.0;
    return
end

% P1^1
P1 = -sqrt(1 - x*x);

if n == 1
    Pn1 = P1;
    return
end

% P2^1
P2 = 3*x*P1;

if n == 2
    Pn1 = P2;
    return
end

Pnm2 = P1;
Pnm1 = P2;

for k = 3:n
    Pn1 = ((2*k-1)*x*Pnm1 - k*Pnm2)/(k-1);
    Pnm2 = Pnm1;
    Pnm1 = Pn1;
end
