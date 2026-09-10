function F = Kaula_INC_dot(l, m, p, INC)
%--------------------------------------------------------------------------
%  Computes dF_{lmp}/dINC, the derivative of Kaula's inclination function
%  with respect to inclination.
%--------------------------------------------------------------------------
%
%   F = Kaula_INC_dot(l, m, p, INC)
%
%   ------
%   Inputs
%   ------
%   l, m, p           degree, order and index of Kaula's function
%   INC       (1,1)   mean inclination [rad]
%
%   -------
%   Outputs
%   -------
%   F         (1,1)   value of dF_{lmp}/dINC []
%
%--------------------------------------------------------------------------

% t limits
if p <= (l-m)/2
    t_max = p;
elseif mod(l-m,2)
    t_max = (l-m-1)/2;
else
    t_max = (l-m)/2;
end

%% Computation of the summation with respect to t
F = 0;
for i = 0:t_max  
    S = 0;
    S_p = 0;
    for s = 0:m
        % c limits
        if m+i >= p+s
            c_min = 0;
        else
            c_min = p-i-m+s;
        end
        if m+p+i >= l+s
            c_max = l-m-2*i+s;
        else
            c_max = p-i;
        end
    
        %% Computation of the summation with respect to c
        C = 0;
        for c = c_min:c_max
            if mod(l-m,2)
                ex = c-0.5*(l-m-1);
            else
                ex = c-0.5*(l-m);
            end
            C = C + factorial(l-m-2*i+s)/(factorial(c)*factorial(l-m-2*i+s-c))*factorial(m-s)/(factorial(p-i-c)*factorial(m-s-p+i+c))/(factorial(l-m-2*i))*(-1)^(ex);
        end
        S = S + factorial(m)/(factorial(s)*factorial(m-s))*cos(INC)^s*C;
        S_p = S_p + factorial(m)/(factorial(s)*factorial(m-s))*s*cos(INC)^(s-1)*sin(INC)*C;
    end
    F = F + factorial(2*l-2*i)/(factorial(i)*factorial(l-i)*2^(2*l-2*i))*(sin(INC).^(l-m-2*i)*S_p + (l-m-2*i)*sin(INC).^(l-m-2*i-1).*cos(INC).*S);
end

end