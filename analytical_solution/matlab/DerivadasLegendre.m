function Pnd = DerivadasLegendre(INC, k)
%DERIVADASLEGENDRE d/dI [P_k(cos(INC))], the derivative of the ordinary
%Legendre polynomial of degree k with respect to the inclination INC
%(not with respect to its argument cos(INC)).
%
%   Pnd = DerivadasLegendre(INC, k)
%   INC : inclination [rad]
%   k   : polynomial degree (2 <= k <= 50)

    maxDegree= 50;
    
    x = cos(INC);
    s = sin(INC);
    
    P = zeros(maxDegree + 1, 1);
    Pd = zeros(maxDegree + 1, 1);
    
    % Initial values of the Legendre polynomials Pn(x)
    P(1) = 1;                % P0
    P(2) = x;                % P1
    
    % Initial values of the derivatives Pn'(x) with respect to x
    Pd(1) = 0;               % P0'
    Pd(2) = 1;               % P1'
    
    % Recurrence for the polynomials and their derivatives with respect to x
    for n = 1:(maxDegree - 1)
        % Standard relation: (n+1)P_{n+1} = (2n+1)xP_n - nP_{n-1}
        P(n+2) = ((2*n + 1)*x*P(n+1) - n*P(n)) / (n+1);
        
        % Relation for the derivative: P'_{n+1} = (n+1)P_n + xP'_n
        Pd(n+2) = (n+1)*P(n+1) + x*Pd(n+1);
    end
    
    % The quantity required is d(Pn(cos(I))) / dI
    % By the chain rule: dPn/dI = dPn/dx * dx/dI = Pd * (-sin(I))
    
    % Build the vector Pnd from n=2 to maxDegree
    n_indices = 2:maxDegree;
    Pnd = -s * Pd(n_indices + 1);

    Pnd = Pnd(k-1);
end