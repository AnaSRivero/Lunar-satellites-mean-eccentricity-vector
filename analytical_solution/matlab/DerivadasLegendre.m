function Pnd = DerivadasLegendre(INC, k)
%DERIVADASLEGENDRE d/dI [P_k(cos(INC))], the derivative of the ordinary
%Legendre polynomial of degree k with respect to the inclination INC
%(not with respect to its argument cos(INC)).
%
%   Pnd = DerivadasLegendre(INC, k)
%   INC : inclination [rad]
%   k   : polynomial degree (2 <= k <= 50)

    gradoMax= 50;
    
    x = cos(INC);
    s = sin(INC);
    
    P = zeros(gradoMax + 1, 1);
    Pd = zeros(gradoMax + 1, 1);
    
    % Valores iniciales para los Polinomios de Legendre Pn(x)
    P(1) = 1;                % P0
    P(2) = x;                % P1
    
    % Valores iniciales para las derivadas Pn'(x) respecto a x
    Pd(1) = 0;               % P0'
    Pd(2) = 1;               % P1'
    
    % Recurrencia para Polinomios y sus derivadas respecto a x
    for n = 1:(gradoMax - 1)
        % Relación estándar: (n+1)P_{n+1} = (2n+1)xP_n - nP_{n-1}
        P(n+2) = ((2*n + 1)*x*P(n+1) - n*P(n)) / (n+1);
        
        % Relación para la derivada: P'_{n+1} = (n+1)P_n + xP'_n
        Pd(n+2) = (n+1)*P(n+1) + x*Pd(n+1);
    end
    
    % Lo que tú buscas es d(Pn(cos(I))) / dI
    % Por regla de la cadena: dPn/dI = dPn/dx * dx/dI = Pd * (-sin(I))
    
    % Creamos el vector Pnd desde n=2 hasta gradoMax
    n_indices = 2:gradoMax;
    Pnd = -s * Pd(n_indices + 1);

    Pnd = Pnd(k-1);
end