function v = corner_speed(kappa, Cl, mu, n)                      % CHANGED: n added
    m = 850;  g = 9.81;  rho = 1.2;  A = 1;
    vmax = 100;                       % cap in m/s (360 kph)
    if nargin < 2, Cl = 3.8; end
    if nargin < 3, mu = 2; end
    if nargin < 4, n = 0; end                                    % NEW
    mu_v = mu * ones(size(kappa));                               % NEW: start with constant grip
    for it = 1:10                                                % NEW: loop wraps your 4 lines
        den = m * kappa - 0.5 * mu_v * rho * Cl * A;             % CHANGED: mu -> mu_v
        v = vmax * ones(size(kappa));
        ok = den > 0;
        v(ok) = min(sqrt(mu_v(ok) * m * g ./ den(ok)), vmax);    % CHANGED: mu -> mu_v
        mu_v = tyre_mu(v, mu, n, m, rho, Cl*A, g);               % NEW: grip at the new speed
    end
end