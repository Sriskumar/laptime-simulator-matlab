function a = accel_max(v, mu, kappa, kf, n, Cl)  
% v in m/s, kappa = path curvature in 1/m, returns max acceleration in m/s^2
% kf = engine force scale factor (default 0.85, calibrated on the real top speed)
if nargin < 2, mu = 2.4; end
if nargin < 3, kappa = 0; end
if nargin < 4, kf = 0.85; end
if nargin < 5, n = 0; end 
if nargin < 6, Cl = 3.8; end  
m = 850;  g = 9.81;  rho = 1.2;  Cd = 1.3;  A = 1;  Fr = 500;
persistent vtab Ftab
if isempty(vtab)
    T = readtable("E:\POLIMI ACADEMIC\PORTFOLIO\Projects\LTS\data\engine_force.csv");
    vtab = T.vCar_kph / 3.6;
    Ftab = T.FxTyre_N;
end
% Engine force at the tyre, constant power beyond the end of the table
Feng = interp1(vtab, Ftab, v);
hi = v > vtab(end);
Feng(hi) = Ftab(end) * vtab(end) ./ v(hi);
Feng = Feng * kf;
% Rear grip: rear weight plus the rear share (62%) of the downforce
Fdown = 0.5 * rho * Cl * A * v.^2;
mu_v  = tyre_mu(v, mu, n, m, rho, Cl*A, g); 
Fgrip = mu_v .* (0.60 * m * g + 0.62 * Fdown);
% Friction circle: the rear tyres also carry 60% of the lateral force
Fy_rear = 0.60 * m * v.^2 .* abs(kappa);
Fgrip   = sqrt(max(Fgrip.^2 - Fy_rear.^2, 0));
Fdrive = min(Feng, Fgrip);
Fdrag  = 0.5 * rho * Cd * A * v.^2;
a = (Fdrive - Fdrag - Fr) / m;
end