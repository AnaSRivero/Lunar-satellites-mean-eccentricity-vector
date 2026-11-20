function F = Kaula_INC_dot(l, m, p, INC)
% Stable dF/dI for even l, m=0, p=l/2 (only case used, in compute_Omega_dot):
% F = P_l(0) P_l(cos I)  ->  dF/dI = -P_l(0) sin(I) P_l'(cos I) = P_l(0) * dP_l(cos I)/dI
if mod(l,2)==0 && m==0 && p==l/2
    F = legendre_rec(l,0) * DerivadasLegendre(INC, l);
else
    error('Kaula_INC_dot (stable): case not implemented');
end
end
