function d = decel_max(v, mu, kappa, n, Cl)     
% v in m/s, kappa = path curvature in 1/m, returns max deceleration (positive) in m/s^2
if nargin < 2, mu = 2.4; end
if nargin < 3, kappa = 0; end
if nargin < 4, n = 0; end 
if nargin < 5, Cl = 3.8; end   
m = 850;  g = 9.81;  rho = 1.2; Cd = 1.3;  A = 1;  Fr = 500;

Fdown = 0.5 * rho * Cl * A * v.^2;
Fdrag = 0.5 * rho * Cd * A * v.^2;

mu_v = tyre_mu(v, mu, n, m, rho, Cl*A, g);

% Friction circle: all four tyres share the lateral force
Fy   = m * v.^2 .* abs(kappa);
Fbrk = sqrt(max((mu_v .* (m * g + Fdown)).^2 - Fy.^2, 0));

d = min((Fbrk + Fdrag + Fr) / m, 5 * g);       % brake system cap at 5 g
end