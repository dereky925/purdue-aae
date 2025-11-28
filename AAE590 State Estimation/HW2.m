clear; clc;close all

linewidth = 3;

% Parse ephemeris data
eph_file = 'ephemeris.txt';

raw = readlines(eph_file);
iStart = find(contains(raw, '$$SOE'), 1, 'first');
iEnd = find(contains(raw, '$$EOE'), 1, 'first');

dataLines = raw((iStart+1):(iEnd-1));
N = numel(dataLines);

tJD = nan(N,1); 
X = nan(N,1); 
Y = nan(N,1); 
Z = nan(N,1);


for k = 1:N
    line = dataLines(k);
    if strlength(strtrim(line))==0
        continue;
    end
    parts = split(line, ',');
    tJD(k) = str2double(strtrim(parts{1}));
    X(k)   = str2double(strtrim(parts{3}));
    Y(k)   = str2double(strtrim(parts{4}));
    Z(k)   = str2double(strtrim(parts{5}));
end

t_sec = (tJD - tJD(1)) * 86400;  % seconds from first sample
t1 = t_sec(1);
t2 = t_sec(end);

% Map to tau in [-1,1]
tau = 2*(t_sec - t1)/(t2 - t1) - 1;

% Chebyshev matrix
chebDesign = @(tauVec, n) buildChebMatrix(tauVec, n);  % returns [N x n]


% Fit for several polynomial degrees
degrees = [3 4 5 6 7 20];  % adjust as desired
RMSE = zeros(numel(degrees),3); % columns for x,y,z

% Pre-allocate storage of fits on a fine grid for plotting
tau_line = linspace(-1,1,2000).';
t_line   = ( (tau_line + 1)/2 )*(t2 - t1) + t1;

fits = struct();
for i = 1:numel(degrees)
    n = degrees(i);
    T = chebDesign(tau, n);

    % Solve LS by QR 
    ax = T \ X;
    ay = T \ Y;
    az = T \ Z;

    % Compute in-sample residuals and RMSE
    rx = X - T*ax; ry = Y - T*ay; rz = Z - T*az;
    RMSE(i,1) = sqrt(mean(rx.^2));
    RMSE(i,2) = sqrt(mean(ry.^2));
    RMSE(i,3) = sqrt(mean(rz.^2));

    % Evaluate on fine grid for plotting
    Tline = chebDesign(tau_line, n);
    fits(i).n   = n;
    fits(i).tx  = Tline*ax;
    fits(i).ty  = Tline*ay;
    fits(i).tz  = Tline*az;
end

% Display RMSE
fprintf('Chebyshev LS fit RMSE [km]\n');
fprintf('  n Deg    RMSE_x       RMSE_y       RMSE_z\n');
for i = 1:numel(degrees)
    fprintf('%5d  %11.4g  %11.4g  %11.4g\n', degrees(i), RMSE(i,1), RMSE(i,2), RMSE(i,3));
end

% Plot
figure('Name','Chebyshev fits vs data','Position',[100 100 1200 1000],'Color','white');

subplot(3,1,1);
plot(t_sec/86400, X, '.', 'DisplayName','data'); hold on;
for i = 1:numel(degrees)
    plot(t_line/86400, fits(i).tx, 'DisplayName',sprintf('Deg: %d', degrees(i)),'LineWidth',linewidth);
end
xlabel('Days from t_1'); ylabel('X (km)'); grid minor; legend show;
set(gca, 'FontSize', 16); 

subplot(3,1,2);

plot(t_sec/86400, Y, '.', 'DisplayName','data'); hold on;
for i = 1:numel(degrees)
    plot(t_line/86400, fits(i).ty, 'DisplayName',sprintf('deg %d', degrees(i)),'LineWidth',linewidth);
end
xlabel('Days from t_1'); ylabel('Y (km)'); grid minor; legend show;
set(gca, 'FontSize', 16);

subplot(3,1,3);

plot(t_sec/86400, Z, '.', 'DisplayName','data'); hold on;
for i = 1:numel(degrees)
    plot(t_line/86400, fits(i).tz, 'DisplayName',sprintf('deg %d', degrees(i)),'LineWidth',linewidth);
end
xlabel('Days from t_1'); ylabel('Z (km)'); grid minor; legend show;
set(gca, 'FontSize', 16);

% Residual plots for a representative degree (pick the last / highest by default)
n_show = degrees(end);
T_show = chebDesign(tau, n_show);
ax = T_show \ X; ay = T_show \ Y; az = T_show \ Z;
rx = X - T_show*ax; ry = Y - T_show*ay; rz = Z - T_show*az;

figure('Name',sprintf('Residuals (degree %d)', n_show),'Position',[100 100 1200 1000],'Color','white');
  


subplot(3,1,1); 
hold on;
grid minor;
plot(t_sec/86400, rx,'LineWidth',linewidth); 
xlabel('Days from t_1'); 
ylabel('X residual (km)'); 
set(gca, 'FontSize', 16); 

subplot(3,1,2); 
hold on;
grid minor;
plot(t_sec/86400, ry,'LineWidth',linewidth); 
xlabel('Days from t_1'); 
ylabel('Y residual (km)'); 
set(gca, 'FontSize', 16);

subplot(3,1,3); 
hold on;
grid minor;
plot(t_sec/86400, rz,'LineWidth',linewidth); 
xlabel('Days from t_1'); 
ylabel('Z residual (km)'); 
set(gca, 'FontSize', 16); 


% Find smallest degree that meets accuracy
tol = 5; % [km] pick RMS target
idx = find(max(RMSE,[],2) < tol, 1, 'first');
if ~isempty(idx)
    fprintf('\nSmallest degree with RMSE < %.3g km on all axes: n = %d\n', tol, degrees(idx));
else
    fprintf('\nNo degree in %s met RMSE < %.3g km on all axes.\n', mat2str(degrees), tol);
end


function T = buildChebMatrix(tauVec, n)

    tauVec = tauVec(:);
    N = numel(tauVec);
    if n < 1
        error('n must be >= 1');
    end
    T = zeros(N, n);
    T(:,1) = 1.0;            % T0
    if n >= 2
        T(:,2) = tauVec;     % T1
    end
    for k = 3:n
        T(:,k) = 2*tauVec.*T(:,k-1) - T(:,k-2);
    end
end



%% Problem 4

clear; clc; close all

LW   = 3;                      % plot line width
FIGP = [100 100 1200 1000];    % figure position
rng(1);                        % noise seed

% Physical/orbital parameters
mu   = 3.986004418e14;         % [m^3/s^2] Earth GM
Re   = 6378e3;                 % [m] Earth mean radius
h    = 450e3;                  % [m] circular orbit altitude
r    = Re + h;                 % [m] orbit radius
nc   = sqrt(mu/r^3);           % [rad/s] mean motion
Torbit = 2*pi/nc;              % [s] orbital period

% Discretization / measurement cadence
dt   = 30;                     % [s] measurement interval
N    = round(Torbit/dt);       % number of steps to ~1 orbit
tf   = N*dt;                   % [s] final time used
t    = (0:N)'*dt;              % time vector, N+1 points

% True initial condition (part b, given)
x0_true = [50; 375; 0.003; -0.01];  % [m, m, m/s, m/s]

% Measurement model
H    = [1 0 0 0;
        0 1 0 0];              % observe planar position only
sigma_meas = 4;                % [m] per-axis 1-sigma
R    = diag([sigma_meas^2, sigma_meas^2]);

Phi_dt = @(dT) cwPhi(dT, nc);  % 4x4 state transition for step dT

% Precompute constant step transition for uniform dt
Phi = Phi_dt(dt);

% Simulate true motion 
x_true = zeros(4, N+1);
x_true(:,1) = x0_true;
for k = 2:(N+1)
    x_true(:,k) = Phi * x_true(:,k-1);
end

% Plot true trajectory
figure('Position',FIGP,'Color','white'); hold on; grid minor; box on;
plot(x_true(1,:), x_true(2,:), '-', 'LineWidth', LW);
xlabel('x (m)'); ylabel('y (m)'); title('True Relative Trajectory');
set(gca, 'FontSize', 16); 
axis equal

% Generate noisy measurements 
z = zeros(2, N+1);
for k = 1:(N+1)
    z(:,k) = H * x_true(:,k) + sigma_meas * randn(2,1);
end

% Plot measurements vs truth
figure('Position',FIGP,'Color','white'); 
tmin = t/60; % minutes
subplot(2,1,1); hold on; grid on; box on;
plot(tmin, x_true(1,:), '-', 'LineWidth', LW, 'DisplayName','x true');
plot(tmin, z(1,:), '.', 'MarkerSize', 10, 'DisplayName','x meas');
xlabel('Time (min)'); ylabel('x (m)'); title('Position Measurements'); legend('Location','best');
set(gca, 'FontSize', 16); 

subplot(2,1,2); hold on; grid on; box on;
plot(tmin, x_true(2,:), '-', 'LineWidth', LW, 'DisplayName','y true');
plot(tmin, z(2,:), '.', 'MarkerSize', 10, 'DisplayName','y meas');
xlabel('Time (min)'); ylabel('y (m)'); legend('Location','best');
set(gca, 'FontSize', 16); 

% Observation matrix
Hk = H;  

% LUMVE from z0 and z1 
A = [H;
     H*Phi];
zb = [z(:,1);
      z(:,2)];
Rb = blkdiag(R, R);

% Weighted least squares
Ri = inv(R);
Wi = blkdiag(Ri, Ri);
P0_hat = inv(A' * Wi * A);
x0_hat = P0_hat * (A' * Wi * zb);

% Report estimation error and covariance
err0 = x0_hat - x0_true;
fprintf('\n(d) LUMVE using z0,z1:\n');
fprintf('x0_hat = [% .6f % .6f % .6f % .6f]^T\n', x0_hat);
fprintf('err0   = [% .6f % .6f % .6f % .6f]^T  |err0|_2 = %.6f\n', err0, norm(err0));
disp('P0_hat (covariance) ='); disp(P0_hat);

% Propagate LUMVE to final time 
% no process noise
Phi_f = cwPhi(tf, nc);              % transition from t0 -> tf
xf_pred = Phi_f * x0_hat;
Pf_pred = Phi_f * P0_hat * Phi_f.'; % Q=0

err_f_pred = xf_pred - x_true(:,end);
fprintf('\n(e) Predicted state at tf from LUMVE (no further measurements):\n');
fprintf('xf_pred = [% .6f % .6f % .6f % .6f]^T\n', xf_pred);
fprintf('err_f_pred = [% .6f % .6f % .6f % .6f]^T  |err|_2 = %.6f\n', err_f_pred, norm(err_f_pred));
disp('Pf_pred ='); disp(Pf_pred);

% Recursive LS (Kalman with Q=0) 
% Start from information provided by (d) but avoid double-counting z1
% Use x0_hat,P0_hat, propagate to k=1 as PRIOR, then begin updates at k=2.
I4 = eye(4);
xhat = zeros(4, N+1);
Phat = zeros(4,4,N+1);

% Prior at k=1
xhat(:,1) = x0_hat;
Phat(:,:,1) = P0_hat;

% Predict to k=1 and optionally update with z1 to avoid double-use of z1,
% we treat k=1 as already "ingested" via LUMVE and START UPDATING AT k=2.
xpred = Phi * xhat(:,1);
Ppred = Phi * Phat(:,:,1) * Phi.';

xhat(:,2)   = xpred;           % store prediction at k=1->2 step
Phat(:,:,2) = Ppred;

for k = 3:(N+1)  % begin at k=3 (t=2*dt). Measurements z(:,k) available every step
    % Predict from k-1 to k
    xpred = Phi * xhat(:,k-1);
    Ppred = Phi * Phat(:,:,k-1) * Phi.';   % Q = 0

    % Kalman update (Recursive LS)
    S = H * Ppred * H.' + R;
    K = Ppred * H.' / S;
    innov = z(:,k) - H * xpred;

    xhat(:,k)   = xpred + K * innov;
    Phat(:,:,k) = (I4 - K*H) * Ppred;    
end

xf_hat = xhat(:,end);
Pf_hat = Phat(:,:,end);
err_f_hat = xf_hat - x_true(:,end);

fprintf('\n(f) Final estimate after Recursive LS (using z_2..z_N):\n');
fprintf('xf_hat = [% .6f % .6f % .6f % .6f]^T\n', xf_hat);
fprintf('err_f_hat = [% .6f % .6f % .6f % .6f]^T  |err|_2 = %.6f\n', err_f_hat, norm(err_f_hat));
disp('Pf_hat ='); disp(Pf_hat);

% Plot
figure('Position',FIGP,'Color','white'); 
subplot(2,2,1); hold on; grid on; box on;
plot(tmin, x_true(1,:), 'k-', 'LineWidth', LW, 'DisplayName','x true');
plot(tmin, xhat(1,:),   'r--', 'LineWidth', LW, 'DisplayName','x est');
xlabel('Time (min)'); ylabel('x (m)');  legend('Location','best');
set(gca, 'FontSize', 16); 

subplot(2,2,2); hold on; grid on; box on;
plot(tmin, x_true(2,:), 'k-', 'LineWidth', LW, 'DisplayName','y true');
plot(tmin, xhat(2,:),   'r--', 'LineWidth', LW, 'DisplayName','y est');
xlabel('Time (min)'); ylabel('y (m)'); legend('Location','best');
set(gca, 'FontSize', 16); 

subplot(2,2,3); hold on; grid on; box on;
plot(tmin, x_true(3,:), 'k-', 'LineWidth', LW, 'DisplayName','xdot true');
plot(tmin, xhat(3,:),   'r--', 'LineWidth', LW, 'DisplayName','xdot est');
xlabel('Time (min)'); ylabel('xdot (m/s)'); legend('Location','best');
set(gca, 'FontSize', 16); 

subplot(2,2,4); hold on; grid on; box on;
plot(tmin, x_true(4,:), 'k-', 'LineWidth', LW, 'DisplayName','ydot true');
plot(tmin, xhat(4,:),   'r--', 'LineWidth', LW, 'DisplayName','ydot est');
xlabel('Time (min)'); ylabel('ydot (m/s)'); legend('Location','best');
sgtitle('State Tracking (Recursive LS)', 'FontWeight','bold');
set(gca, 'FontSize', 16); 

% 2D trajectory: truth vs final estimate
figure('Position',FIGP,'Color','white'); hold on; grid on; box on;
plot(x_true(1,:), x_true(2,:), 'k-', 'LineWidth', LW, 'DisplayName','truth');
plot(xhat(1,:),   xhat(2,:),   'r--', 'LineWidth', LW, 'DisplayName','estimate');
xlabel('x (m)'); ylabel('y (m)'); title('Relative Trajectory: Truth vs Estimate'); legend('Location','best');
set(gca, 'FontSize', 16); 


% CWH equations
function Phi = cwPhi(dT, n)
% cwPhi  State transition for planar CW/Hill equations over interval dT
% State x = [x; y; xdot; ydot], mean motion n [rad/s], psi = n*dT
    psi = n * dT;
    c = cos(psi);
    s = sin(psi);

    % Given blocks
    Phi_rr = [ 4 - 3*c,         0;
               6*(s - psi),     1 ];

    Phi_rv = [ (1/n)*s,               (2/n)*(1 - c);
               (2/n)*(c - 1),         (4/n)*s - (3/n)*psi ];

    Phi_vr = [ 3*n*s,          0;
               6*n*(c - 1),    0 ];

    Phi_vv = [ c,            2*s;
              -2*s,     -3 + 4*c ];

    % Assemble 4x4
    Phi = [Phi_rr, Phi_rv;
           Phi_vr, Phi_vv];
end