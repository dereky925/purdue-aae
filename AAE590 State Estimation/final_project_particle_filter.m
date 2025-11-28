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
% The script is fully self-contained: it synthesizes a DEM and truth
% trajectory, generates noisy LiDAR measurements, runs the PF, and plots
% trajectories and errors so anyone can follow the workflow.

clear; clc; close all;
rng(7); % deterministic seed for repeatability

%% Configuration parameters
cfg.dt          = 0.5;    % [s] time step
cfg.T           = 180;    % [s] total duration
cfg.N           = round(cfg.T / cfg.dt) + 1; % number of epochs
cfg.V           = 55;     % [m/s] nominal groundspeed (constant magnitude)
cfg.omega_cmd   = @(t) 0.01 * sin(0.03*t); % [rad/s] commanded turn rate profile
cfg.h_msl       = @(t) 500 + 20*sin(0.01*t); % [m] UAV altitude (mean sea level)

% DEM (synthetic Grand-Canyon-esque terrain)
cfg.dem.x_extent = 2000; % [m] half-width in x
cfg.dem.y_extent = 2000; % [m] half-width in y
cfg.dem.n_pts    = 250;  % grid resolution (n x n)

% LiDAR / measurement model
cfg.meas.bias     = 1.5;  % [m] constant bias
cfg.meas.sigma    = 3.0;  % [m] 1-sigma noise (Gaussian)
cfg.meas.drop_p   = 0.0;  % probability of dropout (set >0 to test robustness)

% Particle filter settings
cfg.pf.Np          = 1500;   % number of particles
cfg.pf.sigma_vel   = 1.5;    % [m/s] process noise on speed magnitude
cfg.pf.sigma_omega = 0.01;   % [rad/s] process noise on turn rate
cfg.pf.sigma_meas  = cfg.meas.sigma; % assumed measurement noise
cfg.pf.resample_ratio = 0.5; % resample when N_eff < ratio * Np

%% Build synthetic DEM (digital elevation model)
% Elevation model f(x, y) uses additive sinusoids and a canyon trench.
xg = linspace(-cfg.dem.x_extent, cfg.dem.x_extent, cfg.dem.n_pts);
yg = linspace(-cfg.dem.y_extent, cfg.dem.y_extent, cfg.dem.n_pts);
[Xg, Yg] = meshgrid(xg, yg);

% Canyon: negative bump along y-axis plus rolling hills
canyon = -120 * exp(-((Xg/350).^2)) .* exp(-((Yg/900).^2));
hills  = 80 * sin(Xg/500) .* cos(Yg/400) + 40 * sin((Xg+Yg)/450);
H = canyon + hills;

% griddedInterpolant expects NDGRID ordering; transpose meshgrid arrays.
% This matches p(x,y) queries to the same terrain used for plotting.
dem_interp = griddedInterpolant(Xg', Yg', H', "linear", "nearest");

%% Simulate truth trajectory (discrete-time kinematics)
truth.x   = zeros(cfg.N,1);
truth.y   = zeros(cfg.N,1);
truth.psi = zeros(cfg.N,1); % heading [rad], north-east frame
truth.h   = zeros(cfg.N,1); % altitude MSL [m]

for k = 2:cfg.N
    t = (k-1) * cfg.dt;
    omega = cfg.omega_cmd(t); % commanded yaw rate

    % Propagate truth (no process noise in the "truth" model).
    truth.psi(k) = truth.psi(k-1) + omega * cfg.dt;
    truth.x(k)   = truth.x(k-1) + cfg.V * cos(truth.psi(k-1)) * cfg.dt;
    truth.y(k)   = truth.y(k-1) + cfg.V * sin(truth.psi(k-1)) * cfg.dt;
    truth.h(k)   = cfg.h_msl(t);
end

% Ground elevation under truth and LiDAR measurements
truth.ground = dem_interp(truth.x, truth.y);

% LiDAR provides terrain elevation under vehicle (MSL) with bias + noise.
meas = truth.ground + cfg.meas.bias + cfg.meas.sigma * randn(cfg.N,1);
drop_mask = rand(cfg.N,1) < cfg.meas.drop_p;
meas(drop_mask) = NaN; % simulate occasional dropout

%% Initialize particle filter
% State vector per particle: [x; y; psi]. Altitude is assumed known from baro/GNSS.
particles.x   = (rand(cfg.pf.Np,1) - 0.5) * 2 * cfg.dem.x_extent;
particles.y   = (rand(cfg.pf.Np,1) - 0.5) * 2 * cfg.dem.y_extent;
particles.psi = (rand(cfg.pf.Np,1) - 0.5) * 2 * pi;
particles.w   = ones(cfg.pf.Np,1) / cfg.pf.Np; % uniform initial weights

% Storage for estimates
est.x   = zeros(cfg.N,1);
est.y   = zeros(cfg.N,1);
est.psi = zeros(cfg.N,1);

%% Helper for angle-mean to avoid wrap issues
angle_mean = @(theta, w) atan2(sum(w .* sin(theta)), sum(w .* cos(theta)));

%% Particle filter recursion
for k = 2:cfg.N
    t = (k-1) * cfg.dt;
    omega_cmd_k = cfg.omega_cmd(t);

    % 1) Propagate each particle with control + process noise
    dv     = cfg.pf.sigma_vel * randn(cfg.pf.Np,1);
    domega = cfg.pf.sigma_omega * randn(cfg.pf.Np,1);
    v_samp = cfg.V + dv;
    omega_samp = omega_cmd_k + domega;

    particles.psi = particles.psi + omega_samp * cfg.dt;
    particles.x   = particles.x + v_samp .* cos(particles.psi) * cfg.dt;
    particles.y   = particles.y + v_samp .* sin(particles.psi) * cfg.dt;

    % 2) Predict measurement for each particle from DEM
    h_pred = dem_interp(particles.x, particles.y); % terrain elevation MSL

    % 3) Compute importance weights (Eq. 1) using Gaussian likelihood
    z_k = meas(k);
    if ~isnan(z_k)
        residual = z_k - h_pred;
        particles.w = particles.w .* exp(-0.5 * (residual ./ cfg.pf.sigma_meas).^2);
    else
        % If measurement is missing, only carry over prior weights
        particles.w = particles.w;
    end

    % 4) Normalize weights (Eq. 2); guard against degeneracy
    w_sum = sum(particles.w);
    if w_sum <= eps
        % All particles unlikely: re-spawn uniformly and reset weights
        particles.x   = (rand(cfg.pf.Np,1) - 0.5) * 2 * cfg.dem.x_extent;
        particles.y   = (rand(cfg.pf.Np,1) - 0.5) * 2 * cfg.dem.y_extent;
        particles.psi = (rand(cfg.pf.Np,1) - 0.5) * 2 * pi;
        particles.w   = ones(cfg.pf.Np,1) / cfg.pf.Np;
    else
        particles.w = particles.w / w_sum;
    end

    % 5) Compute effective sample size and resample if needed (Eq. 3, Eq. 4)
    Neff = 1 / sum(particles.w.^2);
    if Neff < cfg.pf.resample_ratio * cfg.pf.Np
        particles = systematic_resample(particles);
    end

    % 6) State estimate as weighted mean
    est.x(k)   = sum(particles.x .* particles.w);
    est.y(k)   = sum(particles.y .* particles.w);
    est.psi(k) = angle_mean(particles.psi, particles.w);
end

%% Performance metrics
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

figure;
plot((0:cfg.N-1)*cfg.dt, pos_err, 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Position error [m]');
grid on; title('Position Error vs Time');

figure;
subplot(2,1,1);
plot((0:cfg.N-1)*cfg.dt, truth.ground, 'k', 'LineWidth', 1.2, 'DisplayName','Truth ground');
hold on;
plot((0:cfg.N-1)*cfg.dt, meas, 'c.', 'DisplayName','LiDAR meas');
plot((0:cfg.N-1)*cfg.dt, dem_interp(est.x, est.y), 'r', 'LineWidth', 1.0, 'DisplayName','PF ground est');
ylabel('Terrain elev [m MSL]');
legend; grid on; title('Measurement Track');
subplot(2,1,2);
plot((0:cfg.N-1)*cfg.dt, rad2deg(truth.psi), 'k', 'LineWidth', 1.2, 'DisplayName','Truth heading');
hold on;
plot((0:cfg.N-1)*cfg.dt, rad2deg(est.psi), 'r--', 'LineWidth', 1.0, 'DisplayName','PF heading');
ylabel('Heading [deg]'); xlabel('Time [s]');
legend; grid on;

%% Systematic resampling (low-variance sampler)
function p_out = systematic_resample(p_in)
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
