function C_comp = compute_Omega_dot(a_m, e_m, i_m, param, deg)
% COMPUTE_OMEGA_DOT  Tasa secular de precesión del RAAN (dΩ/dt)
%
%   C_comp = compute_Omega_dot(a_m, e_m, i_m, param, deg)
%
%   Calcula la variación secular del nodo ascendente debida a los armónicos
%   zonales pares del campo gravitatorio lunar (hasta el grado 'deg'),
%   usando elementos medios de Kozai. Los coeficientes Jn se obtienen
%   internamente de CoefMoon.m.
%
%   La evolución del RAAN se obtiene como:
%       Omega(t) = Omega_0 + C_comp * (t - t0)
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
%   a_m     [km]       Semieje mayor medio (escalar)
%   e_m     [-]        Excentricidad media (escalar)
%   i_m     [rad]      Inclinación media (escalar)
%   param   struct     Parámetros del cuerpo central:
%                        param.mu  [km^3/s^2]  Parámetro gravitacional
%                        param.Re  [km]        Radio ecuatorial
%   deg     [-]        Grado máximo de armónicos zonales pares a incluir
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
%   C_comp  [rad/s]    Tasa secular dΩ/dt evaluada en los elementos medios
%
% -------------------------------------------------------------------------
% REFERENCIAS
% -------------------------------------------------------------------------
%   - Kozai, Y. (1959). "The motion of a close earth satellite."
%     AJ, 64, 367.
%   - Brouwer, D. (1959). "Solution of the problem of artificial satellite
%     theory without drag." AJ, 64, 378.
%   - Liu, J.J.F. (1974). "Satellite motion about an oblate Earth."
%     AIAA Journal, 12(11), 1511-1516. DOI: 10.2514/3.49537
%   - Iorio, L. (2003). "The impact of the static part of the Earth's
%     gravity field on some tests of GR with SLR." arXiv:gr-qc/0203050
%
% -------------------------------------------------------------------------
% NOTAS
% -------------------------------------------------------------------------
%   * Todos los inputs deben estar en unidades consistentes (SI o km+s).
%   * Los Jn deben ser los coeficientes del cuerpo central correspondiente
%     (usar GRGM1200A o LP150Q para la Luna, EGM2008 para la Tierra).
%   * La fórmula es singular en la inclinación crítica i ~ 63.4° y sus
%     armónicos superiores (singularidad de Kozai). Evitar su uso en la
%     banda |i - i_crit| < 1°.
%   * J5, J7 (armónicos impares) no contribuyen a la precesión secular
%     del RAAN y no se incluyen.

% =========================================================
% Variables auxiliares
% =========================================================
n = sqrt(1 ./ (a_m/param.Re).^3);   % unperturbed mean motion

% =========================================================
% Bloque J2 (1er orden) + J2^2 (2do orden) + J4 (1er orden)
% + J2^3 y J2*J4 (3er orden) + J3^2/J2
% Kozai (1959), Liu (1974)
% =========================================================
[~,~,~,C_mtx,~] = CoefMoon;
C0 = [0;C_mtx(:,1)];

sum_even = 0;

for i = 2:2:deg
    s = i/2;
    F_dot = Kaula_INC_dot(2*s, 0,s, i_m);
    sum_even = sum_even + F_dot*C0(i)*(1 + i*(i+1)/4*e_m^2);
end

C_comp = n/sin(i_m)/sqrt(1-e_m^2)*sum_even;

end