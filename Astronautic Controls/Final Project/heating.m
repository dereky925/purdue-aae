% Time‐Optimal Powered Descent with Throttle, Fuel, q‐dyn & Chapman Heating
clear; clc; close all
tic

% ==================== Problem Parameters (Earth) ====================
alpha = 0.3;                   % fuel consumption coeff [s/km]
g     = [-9.81e-3; 0; 0];      % gravity in km/s^2 (acting in -r1 direction)
T_min = 1;                     % kg·km/s^2
T_max = 500;                   % kg·km/s^2
gamma_min_deg = 1;             % glide slope angle in degrees

m0  = 2000;                    % initial mass [kg]
z0  = log(m0);                 % log of initial mass

% Discretization parameters
N      = 30;                   % number of intervals
t0     = 0;
r0     = [400;1000;1000];      % initial position [km]
v0     = [0;-5;-5];            % initial velocity [km/s]
gamma_min_rad = gamma_min_deg * pi/180;

% Continuous‐time system matrices
A_ct = zeros(7);
A_ct(1,4)=1; A_ct(2,5)=1; A_ct(3,6)=1;
B_ct = zeros(7,4);
B_ct(4,1)=1; B_ct(5,2)=1; B_ct(6,3)=1;
B_ct(7,4)=-alpha;
c_ct = zeros(7,1);
c_ct(4:6)=g;

% ============== Solver & Bisection Settings =================
ops     = sdpsettings('solver','mosek','verbose',0);
tol     = 1e-3;      % time convergence tolerance [s]
maxIter = 20;        % maximum bisection iterations

% Find feasible upper‐bound T_high by doubling
T_high = N*1.0;  % initial guess
T_low  = 0;
while true
    dt_test = T_high / N;
    [Ad0,Bd0,cd0] = zohOneStep(A_ct, B_ct, c_ct, dt_test);
    [~,~,cons0] = createDecisionVarsAndConstraints(...
        Ad0,Bd0,cd0,dt_test,r0,v0,z0,m0,alpha,T_min,T_max,gamma_min_rad,N);
    sol = optimize(cons0, [], ops);
    if sol.problem==0
        break
    end
    T_high = 2*T_high;
    if T_high > 1e4
        error('Could not find a feasible upper time bound.')
    end
end

% Bisection loop
iterAccumulate = 0;
for iter = 1:maxIter
    iterAccumulate = iterAccumulate + 1;
    T_mid    = (T_low + T_high)/2;
    dt_test  = T_mid / N;
    [Ad1,Bd1,cd1] = zohOneStep(A_ct, B_ct, c_ct, dt_test);
    [~,~,cons1]  = createDecisionVarsAndConstraints(...
        Ad1,Bd1,cd1,dt_test,r0,v0,z0,m0,alpha,T_min,T_max,gamma_min_rad,N);
    sol = optimize(cons1, [], ops);
    if sol.problem==0
        T_high = T_mid;
        bestAd = Ad1;  bestBd = Bd1;   bestcd = cd1;
    else
        T_low = T_mid;
    end
    if (T_high - T_low) < tol
        break
    end
end
fprintf('Total Bisections: %d\n', iterAccumulate);

% Final solve at T_opt
T_opt = T_high;
dt    = T_opt / N;
[X,U,consOpt] = createDecisionVarsAndConstraints(...
    bestAd,bestBd,bestcd,dt,r0,v0,z0,m0,alpha,T_min,T_max,gamma_min_rad,N);
sol = optimize(consOpt, [], ops);
toc

if sol.problem==0
    % Print summary
    fprintf('\n    Minimum time found: %.4f s\n', T_opt);
    xFinal    = value(X{end});
    finalMass = exp(xFinal(7));
    fprintf('    Final mass:      %.4f kg\n', finalMass);
    fuelConsumed = m0 - finalMass;
    fprintf('    Fuel consumed:   %.4f kg\n', fuelConsumed);

    % Unpack trajectories
    Xval = zeros(7, N+1);
    for k=1:N+1
        Xval(:,k) = value(X{k});
    end
    rVals = Xval(1:3, :);    % [r1; r2; r3] in km
    vVals = Xval(4:6, :);    % [v1; v2; v3] in km/s
    sigmaVals = zeros(1,N);
    for k=1:N
        sigmaVals(k) = value(U{k}(4));
    end

    % Atmosphere table (alt [km], rho [kg/m^3])
    atmosTable = [
        0   1.225;
        10  0.4135;
        20  0.0880;
        30  0.0184;
        40  0.00397;
        50  0.00103;
        60  0.0003097;
        70  0.0000828;
        80  0.0000185
    ];
    h_km  = rVals(1,:);
    rho   = interp1(atmosTable(:,1),atmosTable(:,2),h_km,'linear','extrap');
    v_mag = sqrt(sum(vVals.^2,1));    % km/s
    v_mps = v_mag * 1e3;              % m/s

    % Dynamic pressure
    q_dyn = 0.5 * rho .* v_mps.^2;    % Pa
    peakDynamicPressure = max(q_dyn);
    fprintf('    Peak dynamic pressure: %.2f Pa\n', peakDynamicPressure);

    % Chapman heating
    Rn     = 1;            % nose radius [m]
    C_chap = 1.83e-4;      % Chapman heating constant
    q_chap = C_chap * sqrt(rho./Rn) .* v_mps.^3;  % W/m^2
    peakHeating = max(q_chap);
    fprintf('    Peak Chapman heating:   %.2f W/m^2\n', peakHeating);

    % Plotting
    plotLosslessProperty(U, N, dt);
    plotResults(X, U, N, dt, t0, T_opt, T_max);
else
    error('No feasible solution at T_opt = %.4f s', T_opt);
end

% ==================== Local Functions ====================
function [X,U,constraints] = createDecisionVarsAndConstraints(...
    Ad,Bd,cd,dt,r0,v0,z0,m0,alpha,T_min,T_max,gamma_min_rad,N)

    X = cell(N+1,1);  U = cell(N,1);
    for k=1:N+1, X{k} = sdpvar(7,1,'full'); end
    for k=1:N,   U{k} = sdpvar(4,1,'full');   end

    constraints = X{1} == [r0; v0; z0];

    % dynamics
    for k=1:N
        constraints = [constraints, X{k+1} == Ad*X{k} + Bd*U{k} + cd];
    end

    % thrust & mass linearization
    z_lb = zeros(N+1,1);
    for i=0:N
        mass_lb = m0 - alpha*T_max*(i*dt);
        z_lb(i+1) = log(mass_lb);
    end
    for k=1:N
        zk = X{k}(7);
        ak = U{k}(1:3);
        sk = U{k}(4);
        e_lbk   = exp(-z_lb(k));
        lin_app = e_lbk*(1 - (zk - z_lb(k)));
        constraints = [constraints,...
            sk >= T_min * e_lbk,...
            sk <= T_max * lin_app,...
            norm(ak,2) <= sk ];
    end

    % boundary conditions
    constraints = [constraints, X{end}(1:3)==0, X{end}(4:6)==0];

    % glide-slope
    for k=1:N+1
        rk = X{k}(1:3);
        constraints = [constraints, ...
            rk(1) >= tan(gamma_min_rad)*norm(rk(2:3),2)];
    end
end

function [Ad, Bd, cd] = zohOneStep(A,B,c,dt)
    n = size(A,1);
    M = expm(A*dt);
    if rank(A)==n
        Ainv = inv(A);
        IntA = Ainv*(M - eye(n));
    else
        % numerical integral fallback
        steps = 50; dTau = dt/steps;
        IntA = zeros(n);
        for i=1:steps
            tau = (i-1)*dTau;
            IntA = IntA + expm(A*tau)*dTau;
        end
    end
    Ad = M;
    Bd = IntA * B;
    cd = IntA * c;
end

function plotLosslessProperty(U, N, dt)
    diffVals = zeros(1,N);
    for k = 1:N
        a_val     = value(U{k}(1:3));
        sigma_val = value(U{k}(4));
        diffVals(k) = sigma_val - norm(a_val,2);
    end
    t_ctrl = linspace(0, dt*(N-1), N);
    figure('Color','white','Position',[0,100,1000,600]);
    plot(t_ctrl, diffVals, 'k-o','LineWidth',2);
    grid on; hold on;
    xlabel('Time [s]','FontSize',20);
    ylabel('\sigma_k - ||a_k||_2','FontSize',20);
    title('Lossless Property: \sigma_k - ||a_k||_2 vs. Time','FontSize',20);
    set(gca,'FontSize',20);
end

function plotResults(X, U, N, dt, t0, tF, T_max)
    % gather states & controls
    Xval = zeros(7, N+1);
    for k=1:N+1, Xval(:,k) = value(X{k}); end
    rVals = Xval(1:3, :);  vVals = Xval(4:6, :);
    for k=1:N, aVals(:,k) = (vVals(:,k+1)-vVals(:,k))/dt; end
    for k=1:N, sigmaVals(k)= value(U{k}(4));          end

    % 3D trajectory
    t_state = linspace(t0, tF, N+1);
    figure('Color','white','Position',[500,500,1000,800]);
    plot3(-rVals(2,:), rVals(3,:), rVals(1,:), 'b-o','LineWidth',2);
    hold on;
    impact = rVals(:,end);
    scatter3(-impact(2), impact(3), impact(1), 400, 'py','filled','MarkerEdgeColor','k');
    grid on;
    xlabel('r2 [km]','FontSize',20); ylabel('r3 [km]','FontSize',20); zlabel('r1 [km]','FontSize',20);
    title('3D Position Trajectory','FontSize',20);
    legend('Trajectory','Impact Point','FontSize',16,'Location','best');
    set(gca,'FontSize',20); view(75,20); axis equal

    % 2D top view
    figure('Color','white','Position',[500,500,1000,500]);
    plot3(-rVals(2,:), rVals(3,:), rVals(1,:), 'b-o','LineWidth',2);
    grid on; xlabel('r2 [km]','FontSize',20); ylabel('r3 [km]','FontSize',20);
    set(gca,'FontSize',20); view(90,0);

    % positions vs time
    figure('Color','white','Position',[500,0,1000,500]);
    plot(t_state, rVals(1,:),'r-o', t_state, rVals(2,:),'b-o', t_state, rVals(3,:),'g-o','LineWidth',2);
    grid on; xlabel('Time [s]','FontSize',20); ylabel('Position [km]','FontSize',20);
    legend('r1','r2','r3','FontSize',20); title('Position Components vs. Time','FontSize',20);
    set(gca,'FontSize',20);

    % velocities vs time
    figure('Color','white','Position',[3000,1000,1000,500]);
    plot(t_state, vVals(1,:),'r-o', t_state, vVals(2,:),'b-o', t_state, vVals(3,:),'g-o','LineWidth',2);
    grid on; xlabel('Time [s]','FontSize',20); ylabel('Velocity [km/s]','FontSize',20);
    legend('v1','v2','v3','FontSize',20); title('Velocity Components vs. Time','FontSize',20);
    set(gca,'FontSize',20);

    % throttle vs time (fixed)
    t_ctrl = linspace(t0, tF-dt, N);
    throttleLevel = (sigmaVals / T_max) * 100;   % % corrected
    figure('Color','white','Position',[3000,0,1000,500]);
    plot(t_ctrl, throttleLevel,'m-o','LineWidth',2);
    grid on; xlabel('Time [s]','FontSize',20); ylabel('Throttle Level [%]','FontSize',20);
    title('Throttle Level vs. Time','FontSize',20); set(gca,'FontSize',20);

    % control inputs & magnitude
    figure('Color','white','Position',[3000,1500,1000,500]);
    ctrlMag = sqrt(sum(aVals.^2,1));
    plot(t_ctrl, aVals(1,:),'r-','LineWidth',2); hold on;
    plot(t_ctrl, aVals(2,:),'g-','LineWidth',2);
    plot(t_ctrl, aVals(3,:),'b-','LineWidth',2);
    plot(t_ctrl, ctrlMag,'k--','LineWidth',2);
    grid on; xlabel('Time [s]','FontSize',20); ylabel('Control [km/s^2]','FontSize',20);
    legend('a_1','a_2','a_3','mag','FontSize',16,'Location','best');
    title('Control Inputs and Their Magnitude','FontSize',20); set(gca,'FontSize',20);
end