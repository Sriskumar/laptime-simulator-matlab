function mu = tyre_mu(v, mu_ref, n, m, rho, CLA, g)
% Grip coefficient that falls as the load per tyre rises.
% mu_ref : grip at the static load per tyre (m*g/4)
% n      : load sensitivity (n = 0 gives constant grip)
Fz     = (m*g + 0.5*rho*CLA*v.^2) / 4;   % load per tyre at speed v (N)
Fz_ref = m*g / 4;                        % static load per tyre (N)
mu     = mu_ref * (Fz ./ Fz_ref).^(-n);
end