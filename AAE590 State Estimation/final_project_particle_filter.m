% AAE590 Final Project – Terrain-Referenced Particle Filter (MATLAB)
% Derek Yu
%
% This script implements the proposal: use a particle filter to estimate
% a simulated UAV position (x, y, heading) over terrain using noisy LiDAR
% elevation measurements and a digital elevation model (DEM).
%
% Key PF equations referenced in comments:
%  - State propagation (importance proposal): x_k^(i) ~ p(x_k | x_{k-1}, u_k)
%  - Importance weights: w_k^(i) ∝ w_{k-1}^(i) * p(z_k | x_k^(i))         (Eq. 1)
%    Here p(z_k | x_k) is Gaussian: exp(-0.5 * (r/sigma_r)^2).
%  - Weight normalization: w_k^(i) = w_k^(i) / sum_j w_k^(j)               (Eq. 2)
%  - Effective sample size: N_eff = 1 / sum_i (w_k^(i))^2                  (Eq. 3)
%  - Systematic resampling when N_eff < N_thresh                           (Eq. 4)
%
% The script loads a DEM from GeoTIFF (see viz_tif.m), synthesizes a truth
% trajectory over that map, generates noisy LiDAR measurements, runs the PF,
% and plots trajectories and errors so anyone can follow the workflow.

clear; clc; close all;
rng(1); % deterministic seed for repeatability

% State being estimated: x (east), y (north), psi (heading). Altitude is assumed
% known from another sensor. Each particle carries one hypothesis of state plus
% an importance weight that scores how well that hypothesis agrees with the
% LiDAR measurement of ground elevation.

%% Configuration parameters
% Timing and kinematics inputs
cfg.dt          = 0.5;    % [s] time step
cfg.T           = 170;    % [s] total duration (covers ~4.6 km at 27 m/s)
cfg.N           = round(cfg.T / cfg.dt) + 1; % number of epochs
cfg.start.x     = 406788; % [m] UTM Easting start
cfg.start.y     = 3992410; % [m] UTM Northing start
cfg.path.x_target = 411428; % [m] desired end x (eastward)
cfg.V           = (cfg.path.x_target - cfg.start.x) / cfg.T; % [m/s] set to hit target x
cfg.omega_cmd   = @(t) s_turn_profile(t, 80, 0.018); % [rad/s] S-turn then straighten
cfg.h_msl       = @(t) 500 + 20*sin(0.01*t); % [m] UAV altitude (mean sea level)

% DEM (real GeoTIFF)
cfg.dem.tifFile      = 'grand_canyon_rectangle_slice.tin.tif';
cfg.dem.sample_step  = 5;    % decimate raster for speed; set to 1 for full res
cfg.dem.point_step   = 500;  % for optional point cloud visualization

% LiDAR / measurement model
cfg.meas.bias     = 1.5;  % [m] constant bias
cfg.meas.sigma    = 3.0;  % [m] 1-sigma noise (Gaussian)
cfg.meas.drop_p   = 0.0;  % probability of dropout (set >0 to test robustness)
cfg.meas.patch_half = 20; % [m] half-width of square footprint (edge-to-edge = 2*half)
cfg.meas.patch_step = 10; % [m] sample spacing inside the footprint grid

% Particle filter settings (fixed particle count)
cfg.pf.Np           = 100000;  % number of particles
cfg.pf.sigma_vel    = 1.5;   % [m/s] process noise on speed magnitude
cfg.pf.sigma_omega  = 0.01;  % [rad/s] process noise on turn rate
cfg.pf.sigma_meas   = cfg.meas.sigma; % assumed measurement noise
cfg.pf.resample_ratio = 0.5; % resample when N_eff < ratio * Np

%% Load DEM from GeoTIFF (downsampled for runtime)
dem = viz_tif(cfg.dem.tifFile, ...
    'sampleStep', cfg.dem.sample_step, ...
    'pointCloudStep', cfg.dem.point_step, ...
    'makePlots', false);

H = dem.Z;
xg = dem.X(1, :);
yg = dem.Y(:, 1);
dem_interp = dem.interp;
cfg.dem.x_range = dem.bounds.xlim;
cfg.dem.y_range = dem.bounds.ylim;
cfg.dem.center = dem.bounds.center;
x_rng = cfg.dem.x_range;
y_rng = cfg.dem.y_range;
x_span = diff(x_rng);
y_span = diff(y_rng);

% Precompute square patch offsets for LiDAR footprint (cartesian, map-aligned)
patch_axis = -cfg.meas.patch_half : cfg.meas.patch_step : cfg.meas.patch_half;
[patch_dx_grid, patch_dy_grid] = meshgrid(patch_axis, patch_axis);
patch_dx = patch_dx_grid(:)'; % 1 x M offsets in x
patch_dy = patch_dy_grid(:)'; % 1 x M offsets in y
patch_count = numel(patch_dx);

% Keep footprints inside DEM: shrink usable area by patch_half margin
safe_x_rng = [x_rng(1)+cfg.meas.patch_half, x_rng(2)-cfg.meas.patch_half];
safe_y_rng = [y_rng(1)+cfg.meas.patch_half, y_rng(2)-cfg.meas.patch_half];

% Snap start/target into safe area if user-specified values are outside
if cfg.start.x < safe_x_rng(1) || cfg.start.x > safe_x_rng(2) || ...
   cfg.start.y < safe_y_rng(1) || cfg.start.y > safe_y_rng(2)
    warning('Start was outside safe DEM margin. Clipping into interior.');
    cfg.start.x = min(max(cfg.start.x, safe_x_rng(1)), safe_x_rng(2));
    cfg.start.y = min(max(cfg.start.y, safe_y_rng(1)), safe_y_rng(2));
end
if cfg.path.x_target < safe_x_rng(1) || cfg.path.x_target > safe_x_rng(2)
    warning('Path x_target outside safe margin. Clipping into interior.');
    cfg.path.x_target = min(max(cfg.path.x_target, safe_x_rng(1)), safe_x_rng(2));
end
cfg.V = (cfg.path.x_target - cfg.start.x) / cfg.T; % recompute speed if target changed

% Warn if requested start is outside DEM so you notice if simulation leaves map
if cfg.start.x < x_rng(1) || cfg.start.x > x_rng(2) || ...
   cfg.start.y < y_rng(1) || cfg.start.y > y_rng(2)
    warning('Start (%.1f, %.1f) is outside DEM bounds x:[%.1f %.1f], y:[%.1f %.1f]', ...
        cfg.start.x, cfg.start.y, x_rng(1), x_rng(2), y_rng(1), y_rng(2));
end

%% Simulate truth trajectory (discrete-time kinematics)
truth.x   = zeros(cfg.N,1);
truth.y   = zeros(cfg.N,1);
truth.psi = zeros(cfg.N,1); % heading [rad], north-east frame
truth.h   = zeros(cfg.N,1); % altitude MSL [m]
truth.x(1) = cfg.start.x;
truth.y(1) = cfg.start.y;
truth.h(1) = cfg.h_msl(0);
truth.psi(1) = 0; % eastbound
truth_out_of_bounds_warned = false;

for k = 2:cfg.N
    t = (k-1) * cfg.dt;
    omega = cfg.omega_cmd(t); % commanded yaw rate (S-turn then straight)

    % Propagate truth (no process noise in the "truth" model).
    truth.psi(k) = truth.psi(k-1) + omega * cfg.dt;
    truth.x(k)   = truth.x(k-1) + cfg.V * cos(truth.psi(k-1)) * cfg.dt;
    truth.y(k)   = truth.y(k-1) + cfg.V * sin(truth.psi(k-1)) * cfg.dt;
    truth.h(k)   = cfg.h_msl(t);
    % Keep truth footprint inside DEM margin; warn once if clamped
    if truth.x(k) < safe_x_rng(1) || truth.x(k) > safe_x_rng(2)
        truth.x(k) = min(max(truth.x(k), safe_x_rng(1)), safe_x_rng(2));
        if ~truth_out_of_bounds_warned
            warning('Truth x clamped into safe DEM bounds at t=%.1f s', t);
            truth_out_of_bounds_warned = true;
        end
    end
    if truth.y(k) < safe_y_rng(1) || truth.y(k) > safe_y_rng(2)
        truth.y(k) = min(max(truth.y(k), safe_y_rng(1)), safe_y_rng(2));
        if ~truth_out_of_bounds_warned
            warning('Truth y clamped into safe DEM bounds at t=%.1f s', t);
            truth_out_of_bounds_warned = true;
        end
    end
end

% Ground elevation under truth and LiDAR measurements
truth.ground = dem_interp(truth.x, truth.y); % center point (for reference)

% Square-footprint LiDAR: sample a grid around the vehicle, add bias/noise.
% Any samples outside the DEM become NaN, reducing valid points for likelihood.
meas_patch = nan(cfg.N, patch_count);     % each row is one footprint (flattened)
truth_patch_mean = nan(cfg.N,1);          % mean ground over footprint
for k = 1:cfg.N
    pt_x = truth.x(k) + patch_dx;
    pt_y = truth.y(k) + patch_dy;
    in_bounds = pt_x >= x_rng(1) & pt_x <= x_rng(2) & pt_y >= y_rng(1) & pt_y <= y_rng(2);
    patch_truth = dem_interp(pt_x, pt_y);
    patch_truth(~in_bounds) = NaN;
    truth_patch_mean(k) = nanmean(patch_truth(:));
    meas_patch(k,:) = patch_truth + cfg.meas.bias + cfg.meas.sigma * randn(1, patch_count);
end
drop_mask = rand(cfg.N,1) < cfg.meas.drop_p;
meas_patch(drop_mask, :) = NaN; % simulate occasional dropout of full scan
meas_mean = nanmean(meas_patch, 2); % collapsed statistic for plotting

%% Initialize particle filter
% State vector per particle: [x; y; psi]. Altitude is assumed known from baro/GNSS.
% Prior: spread uniformly over DEM bounds (reduced margin) with random heading.
Np = cfg.pf.Np;
particles.x   = safe_x_rng(1) + rand(Np,1) * diff(safe_x_rng);
particles.y   = safe_y_rng(1) + rand(Np,1) * diff(safe_y_rng);
particles.psi = (rand(Np,1) - 0.5) * 2 * pi;
particles.w   = ones(Np,1) / Np; % uniform initial weights

% Storage for estimates (filter outputs at each time)
est.x   = zeros(cfg.N,1);
est.y   = zeros(cfg.N,1);
est.psi = zeros(cfg.N,1);

% Helper for circular mean so headings average correctly across wrap at +/-pi
angle_mean = @(theta, w) atan2(sum(w .* sin(theta)), sum(w .* cos(theta)));

%% Particle filter recursion
diag.Needs = []; % not used; placeholder to collect any future diagnostic flags
diag.weight_sum = zeros(cfg.N,1);
diag.Neff = zeros(cfg.N,1);
diag.num_valid_pts = zeros(cfg.N,1);
diag.resampled = false(cfg.N,1);
diag.nan_weights = false(cfg.N,1);

for k = 2:cfg.N
    t = (k-1) * cfg.dt;
    omega_cmd_k = cfg.omega_cmd(t);

    % 1) Propagate each particle with control + process noise
    %    This is the importance proposal p(x_k | x_{k-1}, u_k).
    Np_curr = numel(particles.w);
    dv     = cfg.pf.sigma_vel * randn(Np_curr,1);
    domega = cfg.pf.sigma_omega * randn(Np_curr,1);
    v_samp = cfg.V + dv;
    omega_samp = omega_cmd_k + domega;

    particles.psi = particles.psi + omega_samp * cfg.dt;
    particles.x   = particles.x + v_samp .* cos(particles.psi) * cfg.dt;
    particles.y   = particles.y + v_samp .* sin(particles.psi) * cfg.dt;

    % 2) Predict square-footprint measurement for each particle from DEM
    %    (same offsets as truth LiDAR footprint).
    patch_x = particles.x + patch_dx; % Np x M via implicit expansion
    patch_y = particles.y + patch_dy; % Np x M
    pred_patch = dem_interp(patch_x, patch_y); % terrain elevation MSL
    inside_patch = (patch_x >= x_rng(1)) & (patch_x <= x_rng(2)) & ...
                   (patch_y >= y_rng(1)) & (patch_y <= y_rng(2));
    pred_patch(~inside_patch) = NaN; % invalidate samples off the map

    % 3) Compute importance weights (Eq. 1) using Gaussian likelihood over patch
    %    Likelihood ~ exp(-0.5 / sigma^2 * sum(residual^2)).
    z_patch = meas_patch(k, :);
    if all(isnan(z_patch))
        % If measurement is missing, only carry over prior weights
        particles.w = particles.w;
    else
        residual = (z_patch - cfg.meas.bias) - pred_patch; % remove bias
        valid = ~isnan(residual);
        sse = sum((residual.^2) .* valid, 2); % sum squared error over valid points
        num_pts = sum(valid, 2);
        diag.num_valid_pts(k) = mean(num_pts); % average over particles (for sanity)
        % Use mean squared error per patch so likelihood doesn't vanish solely
        % because the patch has many samples; also guard divide-by-zero.
        mse = sse ./ max(num_pts, 1);
        log_likelihood = -0.5 * mse ./ (cfg.pf.sigma_meas^2);
        log_likelihood = max(log_likelihood, -50); % floor to prevent underflow
        likelihood = exp(log_likelihood);
        likelihood(num_pts == 0) = 0; % no valid comparisons -> zero likelihood
        particles.w = particles.w .* likelihood;
    end

    % 4) Normalize weights (Eq. 2); guard against degeneracy
    w_sum = sum(particles.w);
    diag.weight_sum(k) = w_sum;
    if w_sum <= eps
        % All particles unlikely: re-spawn uniformly and reset weights
        particles.x   = safe_x_rng(1) + rand(Np_curr,1) * diff(safe_x_rng);
        particles.y   = safe_y_rng(1) + rand(Np_curr,1) * diff(safe_y_rng);
        particles.psi = (rand(Np_curr,1) - 0.5) * 2 * pi;
        particles.w   = ones(Np_curr,1) / Np_curr;
    else
        particles.w = particles.w / w_sum;
    end

    % 5) Compute effective sample size and resample if needed (Eq. 3, Eq. 4)
    %    Systematic resampling keeps diversity when too few particles carry
    %    most of the weight (i.e., weight degeneracy).
    Neff = 1 / sum(particles.w.^2);
    diag.Neff(k) = Neff;
    do_resample = Neff < cfg.pf.resample_ratio * numel(particles.w);
    if do_resample
        particles = systematic_resample(particles);
        diag.resampled(k) = true;
    end

    % 6) State estimate as weighted mean
    est.x(k)   = sum(particles.x .* particles.w);
    est.y(k)   = sum(particles.y .* particles.w);
    est.psi(k) = angle_mean(particles.psi, particles.w);
end

% Predicted patch means from PF trajectory (for plotting)
pf_patch_mean = nan(cfg.N,1);
for k = 1:cfg.N
    pf_patch_mean(k) = nanmean(dem_interp(est.x(k) + patch_dx, est.y(k) + patch_dy));
end

% Quick console diagnostics
bad_epochs = find(diag.weight_sum <= eps);
if ~isempty(bad_epochs)
    fprintf('Warning: weight sum collapsed at epochs: %s\n', mat2str(bad_epochs));
end
nan_est_epochs = find(isnan(est.x) | isnan(est.y) | isnan(est.psi));
if ~isempty(nan_est_epochs)
    fprintf('Warning: NaN estimates at epochs: %s\n', mat2str(nan_est_epochs));
end
if all(isnan(diag.Neff))
    fprintf('Neff is NaN throughout (weights likely collapsing early)\n');
else
    fprintf('Neff min/mean/max: %.1f / %.1f / %.1f\n', nanmin(diag.Neff(2:end)), nanmean(diag.Neff(2:end)), nanmax(diag.Neff(2:end)));
end
fprintf('Avg valid points per particle (likelihood calc) last epoch: %.1f of %d\n', diag.num_valid_pts(end), patch_count);
fprintf('Resampling triggered %d times (out of %d epochs)\n', sum(diag.resampled), cfg.N);

%% Performance metrics
% Position error is Euclidean distance between PF estimate and truth at each time.
pos_err = sqrt( (est.x - truth.x).^2 + (est.y - truth.y).^2 );
rmse    = sqrt(mean(pos_err.^2));

fprintf("Final position error: %.1f m | RMSE over trajectory: %.1f m\n", ...
    pos_err(end), rmse);

%% Visualization
figure;
imagesc(xg, yg, H); axis xy equal;
hold on;
colormap('turbo');
colorbar; title('Terrain DEM with Truth and PF Estimate');
plot(truth.x, truth.y, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Truth');
plot(est.x, est.y, 'r--', 'LineWidth', 1.2, 'DisplayName', 'PF Estimate');
scatter(truth.x(1), truth.y(1), 50, 'go', 'filled', 'DisplayName', 'Start');
legend('Location','best');
xlabel('x [m]'); ylabel('y [m]');
% Note: imagesc flips y by default, so we used axis xy to keep north-up.

figure;
plot((0:cfg.N-1)*cfg.dt, pos_err, 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Position error [m]');
grid on; title('Position Error vs Time');

figure;
subplot(2,1,1);
plot((0:cfg.N-1)*cfg.dt, truth_patch_mean, 'k', 'LineWidth', 1.2, 'DisplayName','Truth patch mean');
hold on;
plot((0:cfg.N-1)*cfg.dt, meas_mean, 'c.', 'DisplayName','LiDAR patch mean');
plot((0:cfg.N-1)*cfg.dt, pf_patch_mean, 'r', 'LineWidth', 1.0, 'DisplayName','PF patch mean est');
ylabel('Terrain elev [m MSL]');
legend; grid on; title('Measurement Track');
subplot(2,1,2);
plot((0:cfg.N-1)*cfg.dt, rad2deg(truth.psi), 'k', 'LineWidth', 1.2, 'DisplayName','Truth heading');
hold on;
plot((0:cfg.N-1)*cfg.dt, rad2deg(est.psi), 'r--', 'LineWidth', 1.0, 'DisplayName','PF heading');
ylabel('Heading [deg]'); xlabel('Time [s]');
legend; grid on;

% Point cloud of all LiDAR patch measurements (all time steps)
x_pc = truth.x(:) + patch_dx; % N x M
y_pc = truth.y(:) + patch_dy; % N x M
z_pc = meas_patch;            % N x M
valid_pc = ~isnan(z_pc);
figure;
scatter3(x_pc(valid_pc), y_pc(valid_pc), z_pc(valid_pc), 6, z_pc(valid_pc), 'filled');
colormap('turbo'); cb = colorbar; cb.Label.String = 'Elevation [m MSL]';
xlabel('x [m]'); ylabel('y [m]'); zlabel('Elevation [m MSL]');
title('LiDAR patch point cloud (all epochs)');
axis equal; grid on; view(45,30);

%% Systematic resampling (low-variance sampler)
function p_out = systematic_resample(p_in)
    % Draw evenly spaced samples through the CDF to reduce variance compared
    % to multinomial resampling. Keeps particles with high weights and drops
    % those with tiny weights, then resets weights to uniform.
    Np = numel(p_in.w);
    cdf = cumsum(p_in.w);
    u0 = rand / Np;
    u = u0 + (0:Np-1)'/Np;
    idx = zeros(Np,1);
    j = 1;
    for i = 1:Np
        while u(i) > cdf(j)
            j = j + 1;
        end
        idx(i) = j;
    end

    p_out.x   = p_in.x(idx);
    p_out.y   = p_in.y(idx);
    p_out.psi = p_in.psi(idx);
    p_out.w   = ones(Np,1) / Np; % weights reset after resampling
end

%% S-turn yaw rate profile: sinusoid for a window, then zero
function omega = s_turn_profile(t, active_duration, omega_amp)
    % Models an "S" maneuver: a sinusoidal yaw rate for active_duration
    % seconds, then straight flight (omega = 0) afterward.
    if t <= active_duration
        omega = omega_amp * sin(2*pi*t/active_duration);
    else
        omega = 0;
    end
end
