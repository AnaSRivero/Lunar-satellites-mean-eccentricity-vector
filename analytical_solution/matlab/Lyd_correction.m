function [e_out, i_out, RAAN_out, AOP_out, M_out] = Lyd_correction(e, i, RAAN, M, AOP, e_sp, i_sp, RAAN_sp, M_sp, AOP_sp)
%LYD_CORRECTION Apply first-order additive Lyddane corrections.

e_out = e + e_sp;
i_out = i + i_sp;
RAAN_out = wrapTo2PiLocal(RAAN + RAAN_sp);
AOP_out = wrapTo2PiLocal(AOP + AOP_sp);
M_out = wrapTo2PiLocal(M + M_sp);

end

function ang = wrapTo2PiLocal(ang)
ang = mod(ang, 2*pi);
ang(ang < 0) = ang(ang < 0) + 2*pi;
end
