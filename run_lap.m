function [vsol, t] = run_lap(mu, n, kappa, vm, ds, Cl)
% RUN_LAP  One full lap (forward + backward pass) for a given grip level.
%   mu     reference friction coefficient (mu_ref)
%   n      load-sensitivity exponent
%   kappa  curvature on the grid (1/m)
%   vm     real speed on the grid (m/s), used for the starting speed
%   ds     grid spacing (m)
%   Cl     lift coefficient (default 3.8)
if nargin < 6, Cl = 3.8; end
N     = numel(vm);
vlim  = corner_speed(abs(kappa), Cl, mu, n);
vmin9 = movmin(vlim, 9);
flip  = [0; double(diff(sign(kappa)) ~= 0)];   % 1 where the curvature changes sign
near  = movmax(flip, 9) > 0;                   % chicane rule: within about 20 m
vlim(near) = vmin9(near);

vf = zeros(N,1);  vf(1) = min(vlim(1), vm(1));
for i = 1:N-1
    a = accel_max(vf(i), mu, kappa(i), 0.85, n, Cl);
    vf(i+1) = min(sqrt(vf(i)^2 + 2*a*ds), vlim(i+1));
end

vb = zeros(N,1);  vb(N) = min(vlim(N), vf(N));
for i = N-1:-1:1
    d = decel_max(vb(i+1), mu, kappa(i+1), n, Cl);
    vb(i) = min(sqrt(vb(i+1)^2 + 2*d*ds), vlim(i));
end

vsol = min(vf, vb);
t    = sum(ds ./ ((vsol(1:end-1) + vsol(2:end)) / 2));
end