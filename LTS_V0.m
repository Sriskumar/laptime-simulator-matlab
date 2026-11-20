%% LTS V0 - quasi-static lap-time simulation, Monza 2025 Q (VER)
% Run from the LTS folder. Telemetry is not included in the repo: put your own
% CSV (Distance, x_m, y_m, v_ms, t_s, Throttle, Brake, nGear) in .\data\
clear all; clc; close all

%% 1. Data and 5 m grid
T  = readtable(fullfile("data","monza_2025Q_VER_fastest.csv"));
s  = (0:5:max(T.Distance))';
x  = interp1(T.Distance, T.x_m,  s);
y  = interp1(T.Distance, T.y_m,  s);
v  = interp1(T.Distance, T.v_ms, s);

sm = s(2:end-1);              % curvature lives on the interior points
vm = v(2:end-1);              % real speed on the same points
ds = sm(2) - sm(1);
N  = numel(sm);

%% 2. Curvature: light smoothing in slow corners, heavy elsewhere
k_fine   = curv(x, y, 7);     % keeps tight corners
k_coarse = curv(x, y, 20);    % suppresses noise on fast sections
slow     = vm < 170/3.6;      % only the chicanes are below 170 kph
kappa    = k_coarse;
kappa(slow) = k_fine(slow);

%% 3. Plausibility cap on curvature
% a car cannot follow a curvature that needs more grip than the model car has at the measured speed
mu_cal      = 2.4;
k_drv=1;
n_load = 0;
mu_cap      =2.4;
ay_lim_real = mu_cap * (9.81 + 0.5*1.2*3.8*vm.^2/850);   % lateral limit (m/s^2) at the real speed
k_allowed   = 1.05 * ay_lim_real ./ vm.^2;               % largest followable curvature (5% tolerance)
capped      = abs(kappa) > k_allowed;
kappa(capped) = sign(kappa(capped)) .* k_allowed(capped);
fprintf("capped %d of %d curvature points\n", sum(capped), numel(kappa))

%% 4. Corner speed limit, with the chicane rule

vlim = corner_speed(abs(kappa), 3.8, mu_cal*k_drv, n_load);
vmin9 = movmin(vlim, 9);
flip  = [0; double(diff(sign(kappa)) ~= 0)];   % 1 where the curvature changes sign
near  = movmax(flip, 9) > 0;                   % within about 20 m of a sign change
vlim(near) = vmin9(near);

%% 5. Forward pass (acceleration)
vf = zeros(N,1);  vf(1) = min(vlim(1), vm(1));
for i = 1:N-1
    a = accel_max(vf(i), mu_cal*k_drv, kappa(i), 0.85, n_load);
    vnext = sqrt(vf(i)^2 + 2*a*ds);
    vf(i+1) = min(vnext, vlim(i+1));
end

%% 6. Backward pass (braking)
vb = zeros(N,1);  vb(N) = min(vlim(N), vf(N));
for i = N-1:-1:1
   d = decel_max(vb(i+1), mu_cal*k_drv, kappa(i+1), n_load);
    vprev = sqrt(vb(i+1)^2 + 2*d*ds);          % highest speed that can still slow to vb(i+1)
    vb(i) = min(vprev, vlim(i));
end

%% 7. Result
vsol    = min(vf, vb);
tvec    = ds ./ ((vsol(1:end-1) + vsol(2:end)) / 2);
t_model = sum(tvec)
t_real  = interp1(T.Distance, T.t_s, sm(end)) - interp1(T.Distance, T.t_s, sm(1))

figure; plot(sm, vlim*3.6, sm, vf*3.6, sm, vsol*3.6, sm, vm*3.6)
xlabel("Distance along lap (m)"); ylabel("Speed (kph)")
legend("corner limit", "forward pass", "model", "actual"); grid on

% shared quantities for the diagnostics below
tr   = interp1(T.Distance, T.t_s, sm);     % real time on the grid
tm   = [0; cumsum(tvec)];                  % model time on the grid
gain = (tr - tr(1)) - tm;                  % model ahead of real, in seconds

%% G1. Identification grid: mu_ref and n against the real lap
mu_list = 2.2:0.1:3.2;
n_list  = 0:0.05:0.25;
err_rms = zeros(numel(mu_list), numel(n_list));   % speed error vs real (m/s)
dt_lap  = zeros(numel(mu_list), numel(n_list));   % model lap minus real lap (s)

for ia = 1:numel(mu_list)
    for ib = 1:numel(n_list)
        [vs, t] = run_lap(mu_list(ia), n_list(ib), kappa, vm, ds);
        err_rms(ia,ib) = sqrt(mean((vs - vm).^2));
        dt_lap(ia,ib)  = t - t_real;
    end
end

[~, idx] = min(err_rms(:));
[ia, ib] = ind2sub(size(err_rms), idx);
fprintf("lowest speed error: mu_ref %.2f, n %.2f, rms %.2f m/s, lap error %+.2f s\n", ...
    mu_list(ia), n_list(ib), err_rms(ia,ib), dt_lap(ia,ib))

figure; imagesc(n_list, mu_list, err_rms); axis xy; colorbar; hold on
contour(n_list, mu_list, dt_lap, [0 0], "w", "LineWidth", 2)
xlabel("n (load sensitivity)"); ylabel("mu_{ref}")
title("Colour: RMS speed error (m/s). White line: lap time matches real")
%% G2. Implied grip vs load (upper envelope of the real car)
CLA = 3.8;  m = 850;  g = 9.81;  rho = 1.2;
ay_real = vm.^2 .* abs(kappa);
mu_imp  = ay_real ./ (g + 0.5*rho*CLA*vm.^2/m);   % grip the real car used
use     = abs(kappa) > 0.005 & ~capped;           % real cornering points, not the capped ones

edges_v = [25 40 55 70 85 100];                   % speed bins (m/s)
vc = nan(1, numel(edges_v)-1);  mu_env = vc;  cnt = vc;
for k = 1:numel(vc)
    sel = use & vm >= edges_v(k) & vm < edges_v(k+1);
    cnt(k) = sum(sel);
    if cnt(k) >= 10
        x = sort(mu_imp(sel));
        mu_env(k) = x(ceil(0.95*numel(x)));       % 95th percentile
        vc(k)     = mean(vm(sel));
    end
end
ok    = ~isnan(mu_env);
ratio = (m*g + 0.5*rho*CLA*vc(ok).^2) / (m*g);    % Fz / Fz_ref at each bin speed
p     = polyfit(log(ratio), log(mu_env(ok)), 1);
n_est = -p(1);  mu_ref_est = exp(p(2));
fprintf("points per bin: %s\n", mat2str(cnt))
fprintf("envelope grip per bin: %s\n", mat2str(mu_env(ok), 3))
fprintf("n = %.3f, mu_ref = %.2f\n", n_est, mu_ref_est)

figure; plot(vc(ok)*3.6, mu_env(ok), "o-"); grid on
xlabel("Speed (kph)"); ylabel("Implied grip, 95th percentile")
%% G3. Calibrate mu_ref to the real lap time for each n
n_set    = [0 0.1 0.2];
mu_match = zeros(size(n_set));
for k = 1:numel(n_set)
    mu_match(k) = fzero(@(mu) lap_err(mu, n_set(k), kappa, vm, ds, t_real), [2.0 3.0]);
    [~, t] = run_lap(mu_match(k), n_set(k), kappa, vm, ds);
    fprintf("n %.2f -> mu_ref %.3f, lap %.4f s (real %.4f)\n", n_set(k), mu_match(k), t, t_real)
end
%% G4. Wing sweep on the three calibrated cars
Cl_list = 3.8 * [0.90 0.95 1.00 1.05 1.10];
t_sweep = zeros(numel(n_set), numel(Cl_list));
for k = 1:numel(n_set)
    for j = 1:numel(Cl_list)
        [~, t_sweep(k,j)] = run_lap(mu_match(k), n_set(k), kappa, vm, ds, Cl_list(j));
    end
end
dt_sweep = t_sweep - t_sweep(:,3);          % lap change vs the real-lap baseline (s)

fprintf("Cl change %%:   %s\n", sprintf("%8.0f", (Cl_list/3.8 - 1)*100))
for k = 1:numel(n_set)
    fprintf("n %.2f dLap s: %s\n", n_set(k), sprintf("%8.3f", dt_sweep(k,:)))
end
figure; plot((Cl_list/3.8 - 1)*100, dt_sweep', "o-"); grid on
xlabel("Cl change (%)"); ylabel("Lap time change (s)"); legend("n = " + string(n_set))
%% G5 - export wing sweep to CSV (results only, no raw telemetry)
n_exp  = [0 0.1 0.2];
mu_exp = [2.233 2.337 2.444];           % calibrated mu_ref from G3
cl_pct = [-10 -5 0 5 10];
res = zeros(numel(n_exp), numel(cl_pct));
for k = 1:numel(n_exp)
    [~, t0] = run_lap(mu_exp(k), n_exp(k), kappa, vm, ds, 3.8);
    for j = 1:numel(cl_pct)
        [~, t] = run_lap(mu_exp(k), n_exp(k), kappa, vm, ds, 3.8*(1 + cl_pct(j)/100));
        res(k,j) = t - t0;
    end
end
Tsw = array2table(res, 'VariableNames', {'Cl_m10','Cl_m5','Cl_0','Cl_p5','Cl_p10'});
Tsw.n = n_exp';  Tsw.mu_ref = mu_exp';
Tsw = movevars(Tsw, {'n','mu_ref'}, 'Before', 1);
disp(Tsw)
if ~isfolder('results'), mkdir('results'); end
writetable(Tsw, fullfile('results','results_wing_sweep.csv'));

%% G6 - wing sweep figure for the README
figure; hold on; grid on
for k = 1:height(Tsw)
    plot(cl_pct, Tsw{k, 3:7}, 'o-', 'DisplayName', sprintf('n = %.1f', Tsw.n(k)));
end
xlabel('Change in Cl (%)'); ylabel('Lap-time change (s)');
title('Downforce sensitivity depends on tyre load sensitivity');
legend('Location','northeast');
try, theme(gcf, 'light'); catch, set(gcf, 'Color', 'w'); end
if ~isfolder('figures'), mkdir('figures'); end
exportgraphics(gcf, fullfile('figures','wing_sweep.png'), 'Resolution', 200);

%% D1. Real vs model time along the lap
figure; plot(sm, -gain)
xlabel("Distance along lap (m)"); ylabel("Real minus model time (s)"); grid on

%% D2. Real braking vs model limit
a_real = -(vm(2:end).^2 - vm(1:end-1).^2) / (2*ds) / 9.81;
a_s    = movmean(a_real, 5);               % light smoothing of the noisy derivative
isbrk  = a_s > 0.5;                        % braking points only
vv     = (30:1:98)';                       % m/s
figure
scatter(vm([false; isbrk])*3.6, a_s(isbrk), 8, "filled"); hold on
plot(vv*3.6, decel_max(vv, mu_cal)/9.81, "LineWidth", 2)
xlabel("Speed (kph)"); ylabel("Deceleration (g)")
legend("real braking points", "model limit"); grid on

%% D3. Grip used by the real car (1 = at the model limit)
ay_real = vm.^2 .* abs(kappa);
util    = ay_real ./ ay_lim_real;
figure
subplot(2,1,1)
plot(sm, util); hold on; yline(1,'--')
ylim([0 1.5]); ylabel('grip used'); xlabel('Distance (m)'); grid on
subplot(2,1,2)
corner = abs(kappa) > 0.01;
plot(vm(corner)*3.6, util(corner), '.'); yline(1,'--')
ylim([0 1.5]); xlabel('Real speed (kph)'); ylabel('grip used'); grid on

%% D4. Gap split: straights vs elsewhere
dg       = [0; diff(gain)];                % seconds gained per 5 m step
straight = abs(kappa) < 0.001;
fprintf("gap total:           %.2f s\n", gain(end));
fprintf("gained on straights: %.2f s\n", sum(dg(straight)));
fprintf("gained elsewhere:    %.2f s\n", sum(dg(~straight)));

%% D5. Zone table: where the gap sits
edges = [sm(1) 800 1100 2000 2300 2900 3800 4300 5000 sm(end)];
names = ["Pit straight" "Rettifilo" "Curva Grande" "Roggia" "Lesmos" ...
         "Serraglio" "Ascari" "Back straight" "Parabolica"];
g_zone = diff(interp1(sm, gain, edges));   % time gained inside each zone
for z = 1:numel(names)
    fprintf("%-14s %5.0f - %5.0f m   %+6.2f s\n", names(z), edges(z), edges(z+1), g_zone(z))
end
fprintf("%-14s %18s %+6.2f s   (gain(end) = %.2f)\n", "Total", "", sum(g_zone), gain(end))

%% D6. Zoom on a window: speed, pedals, gear (change win)
win = [880 1100];
thr = interp1(T.Distance, T.Throttle, sm);
brk = interp1(T.Distance, T.Brake,    sm);
gr  = interp1(T.Distance, T.nGear,    sm, 'previous');
figure
a1 = subplot(3,1,1); plot(sm, vm*3.6, sm, vsol*3.6); ylabel('kph'); legend('real','model'); grid on
a2 = subplot(3,1,2); plot(sm, thr, sm, brk*100); ylabel('% / brake'); legend('throttle','brake x100'); grid on
a3 = subplot(3,1,3); plot(sm, gr); ylabel('gear'); xlabel('Distance (m)'); grid on
linkaxes([a1 a2 a3], 'x'); xlim(win)

%% D7. Straight-line check: 0-100 and 0-200 kph
mu_t = 2;      % traction friction for this test (reference: 2.0 with kf = 1)
kf_t = 1;      % engine force scale for this test (lap uses the accel_max default, 0.85)
ds_s = 1;
x_s  = (0:ds_s:500)';
v_s  = zeros(size(x_s));
for n = 1:numel(x_s)-1
    a = accel_max(v_s(n), mu_t, 0, kf_t);
    v_s(n+1) = sqrt(v_s(n)^2 + 2*a*ds_s);
end
dt   = ds_s ./ ((v_s(1:end-1) + v_s(2:end)) / 2);
t_s  = [0; cumsum(dt)];
t100 = interp1(v_s, t_s, 100/3.6);
t200 = interp1(v_s, t_s, 200/3.6);
fprintf("kf %.2f, mu %.2f:  0-100 kph: %.2f s | 0-200 kph: %.2f s\n", kf_t, mu_t, t100, t200)


%% Local functions 
function kappa = curv(x, y, w)
    xs = smoothdata(x, "gaussian", w);  ys = smoothdata(y, "gaussian", w);
    th = unwrap(atan2(diff(ys), diff(xs)));
    kappa = diff(th) ./ hypot(diff(xs(1:end-1)), diff(ys(1:end-1)));
end
function e = lap_err(mu, n, kappa, vm, ds, t_real)
[~, t] = run_lap(mu, n, kappa, vm, ds);
e = t - t_real;
end