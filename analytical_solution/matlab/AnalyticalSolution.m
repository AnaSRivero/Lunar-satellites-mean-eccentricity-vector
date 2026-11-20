function [e_analy, AOP_analy, RAAN_analy, i_analy, a_analy] = AnalyticalSolution (kep, tau, theta0, deg, a_IC)
%ANALYTICALSOLUTION Analytical mean eccentricity-vector solution.
%
%   [e_analy, AOP_analy, RAAN_analy, i_analy, a_analy] = ...
%       AnalyticalSolution(kep, tau, theta0, deg, a_IC)
%
%   Cook's zonal eccentricity-vector framework extended with the
%   odd-degree tesseral/sectorial forcing terms derived in the paper.
%
%   INPUTS
%     kep    = 5xN osculating Keplerian elements [a; e; i; RAAN; AOP]
%              [km, -, rad, rad, rad] of the reference propagation over
%              one revolution centred at the initial epoch, sampled as in
%              propagate_initial_revolution_campaign.m. They are averaged
%              by mean_elements_first_revolution_v2.m to estimate the mean
%              initial elements of the solution.
%     tau    = dimensionless time vector at which to evaluate the solution
%              (tau = t / sqrt(Re^3/mu), with t measured from the epoch)
%     theta0 = Moon's rotation angle at epoch [rad]
%     deg    = maximum degree of the lunar gravity field
%     a_IC   = initial osculating semi-major axis [km], which defines the
%              length of the averaging revolution
%
%   OUTPUTS
%     e_analy, AOP_analy, RAAN_analy, i_analy, a_analy = mean orbital
%     elements as a function of tau

R_moon = 1737.4;
param.Re = R_moon;               % Moon's radius [km]
param.mu = 4.9028000661637961E+03;    % Moon's gravity [km^3/s^2] (SPICE)

[~,~,~,C_mtx,S_mtx] = CoefMoon;

[kepp, ~, ~] = mean_elements_first_revolution_v2(kep, a_IC);
        
a_m = kepp(1,:);
e_m = kepp(2,:);
i_m = kepp(3,:);
RAAN_m = kepp(4,:);
AOP_m = kepp(5,:);

x0 = [e_m*cos(AOP_m); e_m.*sin(AOP_m)];

%% Secular term of RAAN (even zonal harmonics, see compute_Omega_dot.m)
C_comp = compute_Omega_dot(a_m, e_m, i_m, param, 50);

%Zonal
[K_xi, K_eta, ~] = RACook_const (a_m,i_m./pi.*180,deg);

K = sqrt(abs(K_eta.*K_xi));

% Solution
omegaM = 2.661665148922439e-06*(sqrt(param.Re^3/param.mu));
omegaM = omegaM - C_comp;


n = sqrt(1/(a_m(1)/param.Re)^3);
beta_0 = RAAN_m - theta0;

Q1b = zeros(1,length(tau));
Q2b = zeros(1,length(tau));
if K_xi*K_eta < 0
    % Homogenious
    Phi = [cos(K.*(tau-tau(1))), -K./K_eta.*sin(K.*(tau-tau(1)));
           -K./K_xi.*sin(K.*(tau-tau(1))),  cos(K.*(tau-tau(1)))];
    % Complete solution
    for i = 3:2:deg
        for j = 0:i
            [QC_xi, QS_xi, QC_eta, QS_eta] = coeff_odd(a_m, i_m, n, beta_0, C_mtx, S_mtx, i, j);
            Q1b = Q1b + 1./(K.*(K-j.*omegaM).*(K+j.*omegaM)).*((K.*j.*QS_xi.*omegaM.*cos(j.*omegaM.*tau)-K./K_eta.*QC_eta.*K.^2.*cos(j.*omegaM.*tau)-K.*j.*QC_xi.*omegaM.*sin(j.*omegaM.*tau)-K./K_eta.*QS_eta.*K.^2.*sin(j.*omegaM.*tau)) ...
                +(cos(K.*(tau-tau(1))).*(K.*j.*omegaM.*(-QS_xi.*cos(j.*omegaM.*tau(1))+QC_xi.*sin(j.*omegaM.*tau(1)))+K./K_eta.*K.^2.*(QC_eta.*cos(j.*omegaM.*tau(1))+QS_eta.*sin(j.*omegaM.*tau(1)))) ...
                + sin(K.*(tau-tau(1))).*(K.^2.*(QC_xi.*cos(j.*omegaM.*tau(1))+QS_xi.*sin(j.*omegaM.*tau(1)))+K./K_eta.*j.*K.*omegaM.*(-QC_eta.*sin(j.*omegaM.*tau(1))+QS_eta.*cos(j.*omegaM.*tau(1))))));
            Q2b = Q2b + 1./(K.^3-j.^2.*K.*omegaM.^2).*(((K.^2.*(-K./K_xi).*QC_xi + j.*K.*omegaM.*QS_eta).*cos(j.*omegaM.*tau) + (K.^2.*(-K./K_xi).*QS_xi - j.*K.*omegaM.*QC_eta).*sin(j.*omegaM.*tau)) ...
                -(cos(K.*(tau-tau(1))).*(((-K./K_xi).*K.^2.*QC_xi + j.*K.*omegaM.*QS_eta).*cos(j.*omegaM.*tau(1)) + ((-K./K_xi).*K.^2.*QS_xi - j.*K.*omegaM.*QC_eta).*sin(j.*omegaM.*tau(1)))...
                +sin(K.*(tau-tau(1))).*((-K.^2.*QS_eta-(-K./K_xi).*j.*K.*omegaM.*QC_xi).*sin(j.*omegaM.*tau(1))+((-K./K_xi).*j.*K.*omegaM.*QS_xi-K.^2.*QC_eta).*cos(j.*omegaM.*tau(1)))));
        end
    end
else
    % Homogenious
    Phi = [0.5.*(exp(K.*(tau-tau(1))) + exp(-K.*(tau-tau(1)))), 0.5.*K_xi./K.*(exp(K.*(tau-tau(1))) - exp(-K.*(tau-tau(1))));
           0.5.*K./K_xi.*(exp(K.*(tau-tau(1))) - exp(-K.*(tau-tau(1)))), 0.5.*(exp(K.*(tau-tau(1))) + exp(-K.*(tau-tau(1))))];
    % Complete solution
    for i = 3:2:deg
        for j = 0:i
            [QC_xi, QS_xi, QC_eta, QS_eta] = coeff_odd(a_m, i_m, n, beta_0, C_mtx, S_mtx, i, j);
             Q1b = Q1b + 1./(2.*(K.^2 + j.^2.*omegaM.^2)) .* (-sin(j.*omegaM.*tau(1)).*(QC_xi.*j.*omegaM.*(exp(-K.*(tau-tau(1))) + exp(K.*(tau-tau(1)))) + QS_xi.*K.*(exp(-K.*(tau-tau(1))) - exp(K.*(tau-tau(1)))) + K_xi./K.*QC_eta.*j.*omegaM.*(exp(K.*(tau-tau(1))) - exp(-K.*(tau-tau(1)))) - QS_eta.*K_xi.*(exp(-K.*(tau-tau(1))) + exp(K.*(tau-tau(1))))) ...
                 -cos(j.*omegaM.*tau(1)).*(QC_xi.*K.*(exp(-K.*(tau-tau(1))) - exp(K.*(tau-tau(1)))) - QS_xi.*j.*omegaM.*(exp(-K.*(tau-tau(1))) + exp(K.*(tau-tau(1)))) - QC_eta.*K_xi.*(exp(-K.*(tau-tau(1))) + exp(K.*(tau-tau(1)))) + QS_eta.*j.*omegaM.*K_xi./K.*(exp(-K.*(tau-tau(1))) - exp(K.*(tau-tau(1)))))...
                 +sin(j.*tau.*omegaM).*(2.*QC_xi.*j.*omegaM - 2.*QS_eta.*K_xi)...
                 +cos(j.*tau.*omegaM).*(- 2.*QS_xi.*j.*omegaM - 2.*QC_eta.*K_xi));

            Q2b = Q2b + 1./(2.*(K.^2 + j.^2.*omegaM.^2)).*(-sin(j.*omegaM.*tau(1)).*(QC_eta.*j.*omegaM.*(exp(-K.*(tau-tau(1))) + exp(K.*(tau-tau(1)))) + QS_eta.*K.*(exp(-K.*(tau-tau(1))) - exp(K.*(tau-tau(1)))) + K./K_xi.*QC_xi.*j.*omegaM.*(exp(K.*(tau-tau(1))) - exp(-K.*(tau-tau(1)))) - QS_xi.*K.^2./K_xi.*(exp(-K.*(tau-tau(1))) + exp(K.*(tau-tau(1))))) ...
                 -cos(j.*omegaM.*tau(1)).*(QC_eta.*K.*(exp(-K.*(tau-tau(1))) - exp(K.*(tau-tau(1)))) - QS_eta.*j.*omegaM.*(exp(-K.*(tau-tau(1))) + exp(K.*(tau-tau(1)))) - QC_xi.*K.^2./K_xi.*(exp(-K.*(tau-tau(1))) + exp(K.*(tau-tau(1)))) + QS_xi.*j.*omegaM.*K./K_xi.*(exp(-K.*(tau-tau(1))) - exp(K.*(tau-tau(1)))))...
                 +sin(j.*tau.*omegaM).*(2.*QC_eta.*j.*omegaM - 2.*QS_xi.*K.^2./K_xi)...
                 +cos(j.*tau.*omegaM).*(- 2.*QS_eta.*j.*omegaM - 2.*QC_xi.*K.^2./K_xi));
        end
    end
end
Phi1 = reshape(Phi(:,1:length(tau)),2,1,length(tau));
Phi2 = reshape(Phi(:,(length(tau)+1):end),2,1,length(tau));
Phi = [Phi1,Phi2];

x2 = zeros(2,length(tau));
for i = 1:length(tau)
    x2(:,i) = Phi(:,:,i)*x0 + [Q1b(i);Q2b(i)];
end

e_analy = sqrt(x2(1,:).^2 + x2(2,:).^2);
AOP_analy = atan2(x2(2,:),x2(1,:));
RAAN_analy = RAAN_m + C_comp.*tau;
i_analy = i_m;
a_analy = a_m;
end
