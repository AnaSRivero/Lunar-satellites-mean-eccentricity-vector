function [kepp] = mean2osc_KL_moon (kep, inputtype)

%%-------------------------------------------------------------------------
% First-order mapping (based on the theory developed by Kozai
% and Lyddane) between mean (orbit averaged, with short-period motion (and 
% its mean values with respect to M) removed) orbit elements and 
% osculating (instantaneous) orbit elements.
% 
%
% INPUTS:
%     kep         = vector which contains the classical orbit elements
%                     [ a; e; i; RAAN; AOP; M ]
%     inputtype   = string indicating if input is mean or osculating orbit
%                   elements 
%                     'mean' = mean, 'osc' = osculating
% 
% OUTPUTS:
%     kepp        = vector which contains the classical orbit elements
%                     [ a; e; i; RAAN; AOP; M ]
%
%%-------------------------------------------------------------------------

%  Moon constants
R_moon = 1738;
J2 = 2.032214599835629e-04;

a    = kep(1,:)./R_moon;   % semi-major axis [ ]
e    = kep(2,:);            % eccentricity [ ]
i    = kep(3,:);            % inclination [ rad ]
RAAN = kep(4,:);            % right ascension of the ascending node [ rad ]
AOP  = kep(5,:);            % argument of periapsis [ rad ]
M    = kep(6,:);            % mean anomaly at epoch time t0 [ rad ]
nu = zeros(1, length(M));
for ii = 1:length(M)
    nu(ii)       = KepEqn( M(ii), e(ii), 'rad' );       % true anomaly [ rad ]
end

%  Mod out by 2pi
M = mod(M, 2.*pi);


lambda = sqrt(1-e.^2);
kappa = (sin(i)).^2;
r = a.*(1-e.^2)./(1+e.*cos(nu));
s10 = sin(nu);
s20 = sin(2.*nu);
s30 = sin(3.*nu);
s12 = sin(nu + 2.*AOP);
s1_2 = sin(nu - 2.*AOP);
s22 = sin(2.*nu + 2.*AOP);
s32 = sin(3.*nu + 2.*AOP);
s42 = sin(4.*nu + 2.*AOP);
s52 = sin(5.*nu + 2.*AOP);
c12 = cos(nu + 2.*AOP);
c22 = cos(2.*nu + 2.*AOP);
c32 = cos(3.*nu + 2.*AOP);

% Short-periodic terms
a_sp = J2./(2.*a).*((2-3.*kappa).*(a.^3./r.^3-1./lambda.^3)+3.*kappa.*a.^3./r.^3.*c22);
e_sp = J2.*lambda.^2./(4.*e.*a.^2).*((2-3.*kappa).*(a.^3./r.^3-1./lambda.^3) + 3.*kappa.*a.^3./r.^3.*c22) -...
    3.*J2.*kappa./(4.*e.*a.^2.*lambda.^2).*(c22 + e.*c12 + e.*c32./3) - (J2.*kappa.*e.*(2.*lambda+1).*cos(2.*AOP))./(4.*a.^2.*lambda.^2.*(lambda+1).^2);
i_sp = J2./(8.*a.^2.*lambda.^4).*sin(2.*i).*(3.*c22 + 3.*e.*c12 + e.*c32) - J2.*sin(2.*i).*(2.*lambda.^2-lambda-1).*cos(2.*AOP)./(8.*a.^2.*lambda.^4.*(lambda+1));
RAAN_sp = -3.*J2.*cos(i)./(2.*a.^2.*lambda.^4).*(nu-M+e.*s10-0.5.*(s22 + e.*s12 + e./3.*s32)) - J2.*cos(i).*(2.*lambda.^2-lambda-1).*sin(2.*AOP)./(4.*a.^2.*lambda.^4.*(lambda+1));
AOP_sp = 3.*J2./(2.*a.^2.*lambda.^4).*((4-5.*kappa)./2.*(nu-M+e.*s10) + (5.*kappa - 2)./4.*(s22 + e.*s12 + e./3.*s32) + 1./e.*(((2-3.*kappa)./2).*((1-e.^2./4).*s10 + e./2.*s20 + e.^2./12.*s30) -...
    kappa.*(1./4.*(1+5./4.*e.^2).*s12 - e.^2./16.*s1_2 - 7./12.*(1-e.^2./28).*s32 - 3./8.*e.*s42 - e.^2./16.*s52))) - 3.*J2./(2.*a.^2.*lambda.^4).*(kappa./8 + (1 + 2.*lambda).*(2.*kappa.*lambda.^2 - lambda.^2 - kappa + 1)./(6.*(lambda + 1).^2)).*sin(2.*AOP); 
e_mM_sp = -3.*J2./(2.*a.^2.*lambda.^3).*((2 - 3.*kappa)./2.*((1 - e.^2./4).*s10 + e./2.*s20 + e.^2./12..*s30) - kappa.*(1./4.*(1 + 5./4.*e.^2).*s12 - e.^2./16..*s1_2 - ...
    7./12.*(1-e.^2./28).*s32 - 3./8.*e.*s42 - e.^2./16.*s52)) + e.*J2.*kappa.*(4.*lambda.^3 - lambda.^2 - 18.*lambda - 9).*sin(2.*AOP)./(16.*a.^2.*lambda.^3.*(lambda + 1).^2);
M_sp = e_mM_sp./e;


%% Lyddane's correntions
if(strcmp(inputtype, 'mean'))
    a_osc = (a + a_sp)*R_moon;
    [e_osc, i_osc, RAAN_osc, AOP_osc, M_osc] = Lyd_correction  (e, i, RAAN, M, AOP, e_sp, i_sp, RAAN_sp, M_sp, AOP_sp);
else
    a_osc = (a - a_sp)*R_moon;
     [e_osc, i_osc, RAAN_osc, AOP_osc, M_osc] = Lyd_correction  (e, i, RAAN, M, AOP, -e_sp, -i_sp, -RAAN_sp, -M_sp, -AOP_sp);
end

for kk = 1:length(RAAN_osc)
    if RAAN_osc(kk) < 0 && RAAN(kk) > 0
        RAAN_osc(kk) = RAAN_osc(kk) + 2.*pi;
    end
end

%  Transformed orbit element vector
kepp = [ a_osc ; e_osc ; i_osc ; RAAN_osc ; AOP_osc ; M_osc ];

end