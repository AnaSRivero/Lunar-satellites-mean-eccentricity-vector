function F = Kaula_INC(l, m, p, INC)
% Numerically stable F_{lmp}(INC) for the index pairs used by the solution:
%   odd l, p = (l-1)/2 or (l+1)/2  -> Gooding (2008) closed forms f^{+-}_{lm,1}
%   even l, m = 0, p = l/2          -> F = P_l(0) P_l(cos I)
% Falls back to the original series otherwise.
x = cos(INC); S = sin(INC);
if mod(l,2)==1 && (p == (l-1)/2 || p == (l+1)/2)
    B  = (-1)^floor((l-m+1)/2) * exp(gammaln(l) - (l-1)*log(2) - gammaln((l+1)/2) - gammaln((l+3)/2));
    P  = legendre(l, x);              % P_l^m(x), m = 0..l, Condon-Shortley phase
    Pm = P(m+1); if m < l, Pm1 = P(m+2); else, Pm1 = 0; end
    if m == 0
        fp = 0;  fm = B*Pm1;          % f^-_{l0} = -B S P_l' = B P_l^1
    else
        if m == 1
            PmS = -dPl(l, x);          % P_l^1/sin(I) = -P_l'(cos I): exact, also at I = 0, pi
        else
            PmS = Pm/S;                % ~ sin(I)^(m-1) -> 0 at I = 0 for m >= 2
            if S == 0, PmS = 0; end
        end
        fp = m*B*(-1)^m*PmS;
        fm = B*(-1)^m*(m*x*PmS + Pm1);
    end
    if p == (l-1)/2, F = (fp+fm)/2; else, F = (fp-fm)/2; end
elseif mod(l,2)==0 && m==0 && p==l/2
    F = legendre_rec(l,0)*legendre_rec(l,x);
else
    F = Kaula_INC_series(l,m,p,INC);
end
end

function d = dPl(l, x)
% derivative of the Legendre polynomial P_l at x (stable recurrence)
P0 = 1; P1 = x; d0 = 0; d1 = 1;
if l == 0, d = 0; return; elseif l == 1, d = 1; return; end
for n = 1:l-1
    P2 = ((2*n+1)*x*P1 - n*P0)/(n+1);
    d2 = (n+1)*P1 + x*d1;
    P0 = P1; P1 = P2; d0 = d1; d1 = d2;
end
d = d1;
end
