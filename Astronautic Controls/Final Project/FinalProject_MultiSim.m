% untitled8.m
% ==================== Initialization ====================
clear; clc; close all; tic
% make all plot fonts size 20
set(groot,'defaultAxesFontSize',20,'defaultTextFontSize',20);

% ==================== Problem Parameters (Earth Model) ====================
alpha           = 0.0;           % Fuel consumption coefficient [s/km]
g               = [-9.81e-3;     % Gravity vector [km/s^2]
                    0;
                    0];
T_min           = 1;             % Minimum throttle [kg·km/s^2]
T_max           = 500;           % Maximum throttle [kg·km/s^2]
gamma_min_deg   = 1;             % Minimum glide‐slope angle [deg]
m0              = 2000;          % Initial interceptor mass [kg]
z0              = log(m0);       % Log‐mass state

% ==================== Interceptor MPC Discretization ====================
N               = 30;            % Number of control intervals
r0_1            = [400;1000;1000];       % Interceptor 1 init pos [km]
r0_2            = r0_1 + [0;100;250];        % Interceptor 2 offset in r2,r3
v0              = [0;-4;-4];             % Initial velocities [km/s]
gamma_rad       = gamma_min_deg * pi/180;
t_delay         = 20;            % Delay before interceptor launch [s]

% ==================== ICBM Boost‐Phase Simulation ====================
t_boost         = 180;           % Total boost‐phase duration [s]
N_icbm          = 300;           % Number of samples
pitch_max       = 45 * pi/180;   % Max pitch‐over [rad]
a_mag           = 0.02;          % Thrust accel magnitude [km/s^2]
dt_icbm         = t_boost / (N_icbm-1);

t_icbm    = linspace(0, t_boost, N_icbm);
r_icbm1   = zeros(3, N_icbm);
v_icbm1   = zeros(3, N_icbm);
r_icbm2   = zeros(3, N_icbm);
v_icbm2   = zeros(3, N_icbm);
r_icbm1(:,1) = [0;0;0];
r_icbm2(:,1) = [0;30;80];

for i = 2:N_icbm
    tau      = t_icbm(i);
    pitch    = pitch_max * (tau / t_boost);
    a_th     = [a_mag*cos(pitch); a_mag*sin(pitch); 0];
    a_tot    = a_th + g;
    v_icbm1(:,i) = v_icbm1(:,i-1) + a_tot * dt_icbm;
    r_icbm1(:,i) = r_icbm1(:,i-1) + v_icbm1(:,i-1) * dt_icbm;
    v_icbm2(:,i) = v_icbm2(:,i-1) + a_tot * dt_icbm;
    r_icbm2(:,i) = r_icbm2(:,i-1) + v_icbm2(:,i-1) * dt_icbm;
end

% ==================== Interceptor Continuous‐Time Dynamics ====================
A_ct = zeros(7);
A_ct(1,4)=1; A_ct(2,5)=1; A_ct(3,6)=1;
B_ct = zeros(7,4);
B_ct(4,1)=1; B_ct(5,2)=1; B_ct(6,3)=1; B_ct(7,4)=-alpha;
c_ct = zeros(7,1);
c_ct(4:6) = g;  % make c_ct 7×1

% ==================== Solver & Bisection Setup ====================
ops    = sdpsettings('solver','mosek','verbose',0);
tol    = 1e-3;
maxIt  = 20;

T_low   = 0;
T_high  = t_boost - t_delay;   % intercept must occur during boost
bestAd = []; bestBd = []; bestcd = [];

for iter = 1:maxIt
    T_mid = 0.5*(T_low + T_high);
    dt_mid= T_mid / N;
    [Ad,Bd,cd] = zohOneStep(A_ct,B_ct,c_ct,dt_mid);
    t_int      = min(t_boost, t_delay + T_mid);
    r_t1       = interp1(t_icbm', r_icbm1', t_int)';
    r_t2       = interp1(t_icbm', r_icbm2', t_int)';
    [X1,U1,X2,U2,cons] = createDecisionVarsAndConstraintsTwo(...
        Ad,Bd,cd,dt_mid, r0_1,r0_2,v0,z0,m0, ...
        alpha,T_min,T_max,gamma_rad,N, r_t1,r_t2);
    sol = optimize(cons, [], ops);
    if sol.problem==0
        T_high = T_mid;
        bestAd = Ad; bestBd = Bd; bestcd = cd;
    else
        T_low  = T_mid;
    end
    if (T_high - T_low) < tol, break; end
end

% Final solve
T_opt        = T_high;
dt           = T_opt / N;
t_int_global = min(t_boost, t_delay + T_opt);
r_int_pt1    = interp1(t_icbm',r_icbm1',t_int_global)';
r_int_pt2    = interp1(t_icbm',r_icbm2',t_int_global)';

[X1,U1,X2,U2,consOpt] = createDecisionVarsAndConstraintsTwo(...
    bestAd,bestBd,bestcd,dt, r0_1,r0_2,v0,z0,m0, ...
    alpha,T_min,T_max,gamma_rad,N, r_int_pt1,r_int_pt2);
sol = optimize(consOpt, [], ops);
if sol.problem~=0
    error('No feasible solution at T_opt = %.4f s', T_opt);
end

toc

% ==================== Display Key Results ====================
fprintf('Interceptor burn time: %.2f s\n', T_opt);
fprintf('Global intercept time: %.2f s\n\n', t_int_global);
Isp = 300; g0 = 9.81;
s1 = cellfun(@(u)value(u(4)), U1);
s2 = cellfun(@(u)value(u(4)), U2);
fuel1 = sum((s1*1e3)/(Isp*g0)) * dt;
fuel2 = sum((s2*1e3)/(Isp*g0)) * dt;
fprintf('Fuel consumed (Int 1): %.2f kg\n', fuel1);
fprintf('Fuel consumed (Int 2): %.2f kg\n\n', fuel2);
Xv1 = zeros(7,N+1); Xv2 = zeros(7,N+1);
for k=1:N+1, Xv1(:,k)=value(X1{k}); Xv2(:,k)=value(X2{k}); end
vpk1 = max(sqrt(sum(Xv1(4:6,:).^2,1)));
vpk2 = max(sqrt(sum(Xv2(4:6,:).^2,1)));
Tstag1 = 288 + (vpk1*1e3)^2/(2*1004);
Tstag2 = 288 + (vpk2*1e3)^2/(2*1004);
fprintf('Peak T_{stag} (Int 1): %.2f K\n', Tstag1);
fprintf('Peak T_{stag} (Int 2): %.2f K\n\n', Tstag2);

% ==================== Plots ====================
[~, idx_int] = min(abs(t_icbm - t_int_global));

% 1) 3D Trajectories: two interceptors vs. two ICBMs, post‐intercept in gray
figure('Color','white','Position',[100,100,900,700]); hold on; grid on; axis equal;
plot3(-Xv1(2,:),Xv1(3,:),Xv1(1,:),'b-o','LineWidth',2);
plot3(-Xv2(2,:),Xv2(3,:),Xv2(1,:),'g-s','LineWidth',2);
% ICBM before intercept
plot3(-r_icbm1(2,1:idx_int),r_icbm1(3,1:idx_int),r_icbm1(1,1:idx_int),'r-','LineWidth',2);
plot3(-r_icbm2(2,1:idx_int),r_icbm2(3,1:idx_int),r_icbm2(1,1:idx_int),'m-','LineWidth',2);
% ICBM after intercept (gray)
plot3(-r_icbm1(2,idx_int:end),r_icbm1(3,idx_int:end),r_icbm1(1,idx_int:end),...
      'Color',[0.5 0.5 0.5],'LineWidth',2);
plot3(-r_icbm2(2,idx_int:end),r_icbm2(3,idx_int:end),r_icbm2(1,idx_int:end),...
      'Color',[0.5 0.5 0.5],'LineWidth',2);
% intercept points
plot3(-r_int_pt1(2),r_int_pt1(3),r_int_pt1(1),'kp','MarkerFaceColor','y','MarkerSize',12);
plot3(-r_int_pt2(2),r_int_pt2(3),r_int_pt2(1),'kd','MarkerFaceColor','c','MarkerSize',12);

xlabel('r_2 [km]','FontWeight','bold');
ylabel('r_3 [km]','FontWeight','bold');
zlabel('r_1 [km]','FontWeight','bold');
title('3D Trajectories: Two Interceptors vs. Two ICBMs','FontWeight','bold');
legend('Interceptor 1','Interceptor 2', ...
       'ICBM 1 (boost)','ICBM 2 (boost)', ...
       'ICBM 1 (post)','ICBM 2 (post)', ...
       'Intercept Pt 1','Intercept Pt 2', ...
       'Location','best','FontSize',20);
view(45,25); camproj('perspective'); axis tight; rotate3d on;

% 2) Relative distance between interceptors
t_state = linspace(t_delay, t_delay+T_opt, N+1);
dist_rel = sqrt(sum((Xv1(1:3,:) - Xv2(1:3,:)).^2,1));
figure('Color','white','Position',[200,200,800,500]);
plot(t_state, dist_rel,'k-o','LineWidth',2); grid on;
xlabel('Time [s]'); ylabel('Separation [km]');
title('Relative Distance Between Interceptors','FontWeight','bold');

% 3) Lossless property
plotLosslessProperty(U1, N, dt);
plotLosslessProperty(U2, N, dt);

% 4) State & control histories
plotResults(X1, U1, N, dt, t_delay, t_int_global, T_max);
plotResults(X2, U2, N, dt, t_delay, t_int_global, T_max);


% ==================== Local Function Definitions ====================
function [X1,U1,X2,U2,cons] = createDecisionVarsAndConstraintsTwo(...
    Ad,Bd,cd,dt, r0_1,r0_2,v0,z0,m0, ...
    alpha,T_min,T_max,gamma,N, r_t1,r_t2)

    X1 = cell(N+1,1); U1 = cell(N,1);
    X2 = cell(N+1,1); U2 = cell(N,1);
    for k=1:N+1
        X1{k}=sdpvar(7,1,'full');  X2{k}=sdpvar(7,1,'full');
    end
    for k=1:N
        U1{k}=sdpvar(4,1,'full');  U2{k}=sdpvar(4,1,'full');
    end

    cons = [ X1{1}==[r0_1;v0;z0],  X2{1}==[r0_2;v0;z0] ];
    for k=1:N
        cons = [cons,
            X1{k+1}==Ad*X1{k}+Bd*U1{k}+cd,
            X2{k+1}==Ad*X2{k}+Bd*U2{k}+cd ];
    end

    z_lb = arrayfun(@(i) log(m0 - alpha*T_max*(i*dt)), 0:N)';
    for k=1:N
        z1 = X1{k}(7); a1=U1{k}(1:3); s1=U1{k}(4);
        elb=exp(-z_lb(k)); lap=elb*(1-(z1-z_lb(k)));
        cons=[cons, s1>=T_min*elb, s1<=T_max*lap, norm(a1,2)<=s1];

        z2 = X2{k}(7); a2=U2{k}(1:3); s2=U2{k}(4);
        lap=elb*(1-(z2-z_lb(k)));
        cons=[cons, s2>=T_min*elb, s2<=T_max*lap, norm(a2,2)<=s2];
    end

    cons=[cons,
        X1{end}(1:3)==r_t1,  X2{end}(1:3)==r_t2 ];

    for k=1:N+1
        r1k=X1{k}(1:3); r2k=X2{k}(1:3);
        cons=[cons,
            r1k(1)>=tan(gamma)*norm(r1k(2:3),2),
            r2k(1)>=tan(gamma)*norm(r2k(2:3),2)];
    end
end

function [Ad,Bd,cd] = zohOneStep(A,B,c,dt)
    M=expm(A*dt);
    if rank(A)==size(A,1)
        IntA=A\(M-eye(size(A)));
    else
        steps=50; d=dt/steps; IntA=zeros(size(A));
        for i=1:steps
            IntA=IntA+expm(A*(i-1)*d)*d;
        end
    end
    Ad=M; Bd=IntA*B; cd=IntA*c;
end

function plotLosslessProperty(U,N,dt)
    diffVals=zeros(1,N);
    for k=1:N
        diffVals(k)=value(U{k}(4))-norm(value(U{k}(1:3)),2);
    end
    t=linspace(0,dt*(N-1),N);
    figure('Color','white','Position',[0,100,1000,600]);
    plot(t,diffVals,'k-o','LineWidth',2); grid on;
    xlabel('Time [s]'); ylabel('\sigma_k - ||a_k||_2');
    title('Lossless Property vs. Time');
end

function plotResults(X,U,N,dt,t0,tF,T_max)
    Xmat=zeros(7,N+1);
    for k=1:N+1, Xmat(:,k)=value(X{k}); end
    r=Xmat(1:3,:); v=Xmat(4:6,:);
    a=diff(v,1,2)/dt;
    s=cellfun(@(u)value(u(4)),U)/T_max*100;
    t_s=linspace(t0,tF,N+1);
    t_c=linspace(t0,tF-dt,N);

    figure('Color','white','Position',[500,0,1000,500]);
    plot(t_s,r,'-o','LineWidth',2); grid on;
    xlabel('Time [s]'); ylabel('Position [km]');
    legend('r_1','r_2','r_3'); title('Interceptor Positions vs. Time');

    figure('Color','white','Position',[300,200,1000,500]);
    plot(t_s,v,'-o','LineWidth',2); grid on;
    xlabel('Time [s]'); ylabel('Velocity [km/s]');
    legend('v_1','v_2','v_3'); title('Interceptor Velocities vs. Time');

    figure('Color','white','Position',[300,400,1000,500]);
    plot(t_c,s,'m-o','LineWidth',2); grid on;
    xlabel('Time [s]'); ylabel('Throttle [%]');
    title('Interceptor Throttle vs. Time');

    figure('Color','white','Position',[300,600,1000,500]);
    plot(t_c,[a; sqrt(sum(a.^2,1))],'-','LineWidth',2); grid on;
    xlabel('Time [s]'); ylabel('Acceleration [km/s^2]');
    legend('a_1','a_2','a_3','||a||','Location','best');
    title('Interceptor Control Inputs vs. Time');
end