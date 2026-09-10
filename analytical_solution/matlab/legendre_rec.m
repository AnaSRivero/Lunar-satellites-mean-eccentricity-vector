function Pn = legendre_rec(n, x)
%#codegen

Pn = 0.0; % Inicialización obligatoria

P0 = 1.0;
P1 = x;

if n == 0
    Pn = P0;
    return
elseif n == 1
    Pn = P1;
    return
end

for k = 2:n
    Pn = ((2*k-1)*x*P1 - (k-1)*P0)/k;
    P0 = P1;
    P1 = Pn;
end
