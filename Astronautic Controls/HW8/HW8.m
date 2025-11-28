% =========================================================================
% 
% Filename:       HW8.m
% Author:         Derek Yu
% Institution:    Purdue University
% Course:         AAE590 - Applied Control in Astronautics
% Professor:      Dr. Kenshiro Oguri
% Assignment:     HW 8
% Semester:       Spring 2025
% 
% Description: Homework 8
%
% =========================================================================

%% Problem 8

clear; close all; clc;

% Set flag for saturation of control input (true = apply saturation, false = no saturation)
flag_saturate = true;

% Saturation threshold for control input in dimensional units (km/s^2); 1 m/s^2 = 0.001 km/s^2.
u_sat_dim = 0.001;

% Dimensional gravitational parameter (km^3/s^2)
mu_Earth = 3.9860e5;       

% Actual Chief orbital elements (from table)
a_chief_actual     = 10000;    % [km]
e_chief     = 0.3;
i_chief     = 50;            % [deg]
Omega_chief = 50;            % [deg]
omega_chief = 40;            % [deg]
nu_chief    = 0;             % [deg]

% Actual Deputy orbital elements (from table)
a_deputy_actual     = 11500;   % [km]
e_deputy     = 0.15;
i_deputy     = 35;           % [deg]
Omega_deputy = 50;           % [deg]
omega_deputy = 40;           % [deg]
nu_deputy    = 0;            % [deg]

% Define scaling factors
L_scale = a_chief_actual;                  % Length scale in km (set chief semi-major axis to 1)
T_scale = sqrt(a_chief_actual^3 / mu_Earth); % Time scale in seconds
V_scale = L_scale / T_scale;               % Velocity scale in km/s
A_scale = L_scale / T_scale^2;             % Acceleration scale in km/s^2

% Convert simulation time to nondimensional units
T_period_actual = 2*pi * sqrt(a_chief_actual^3 / mu_Earth);  
T_period_nd = T_period_actual / T_scale;     % should be 2*pi
tspan_nd = [0, 5 * T_period_nd];              % [0, 10*pi] for 5 revolutions

% Compute Chief and Deputy States in Dimensional Units
[rC_eci_dim, vC_eci_dim] = oe2rv(a_chief_actual, e_chief, i_chief, Omega_chief, omega_chief, nu_chief, mu_Earth);
[rD_eci_dim, vD_eci_dim] = oe2rv(a_deputy_actual, e_deputy, i_deputy, Omega_deputy, omega_deputy, nu_deputy, mu_Earth);

% Nondimensionalize the Chief and Deputy States:
rC_eci = rC_eci_dim / L_scale;
vC_eci = vC_eci_dim / V_scale;
rD_eci = rD_eci_dim / L_scale;
vD_eci = vD_eci_dim / V_scale;

% Form Initial Relative State in Nondimensional Units:
dr0 = rD_eci - rC_eci;
dv0 = vD_eci - vC_eci;
x0  = [dr0; dv0];

% Define Control Gains (Nondimensional):
Kr = eye(3);  % identity gain
P_scenarios   = {0.2*eye(3), 2*eye(3), 20*eye(3)};
scenarioNames = {'Under-Damped','Critically-Damped','Over-Damped'};

% Define figure size for all figures:
figW = 1100; 
figH = 800;

% Loop Over Damping Cases and Integrate:
for idx = 1:3
    P = P_scenarios{idx};
    scenarioTag = scenarioNames{idx};
    
    % Set x offset based on case:
    if idx == 1
        x_offset = 0;
    elseif idx == 2
        x_offset = 800;
    elseif idx == 3
        x_offset = 2600;
    end
    
    % Define y offsets for the figures:
    y_offsets = [500 400 300 200 100];
    
    % Integrate nondimensional closed-loop relative dynamics:
    options = odeset('AbsTol',1E-12,'RelTol',1E-12);
    [t_nd, xSol_nd] = ode45(@(t,x) closedLoopDynamics(x,Kr,P), tspan_nd, x0, options);
    
    % Extract relative state components (nondimensional):
    dr_nd = xSol_nd(:,1:3);
    dv_nd = xSol_nd(:,4:6);
    
    % Compute control input in nondimensional units:
    uVals = zeros(length(t_nd),3);
    for k = 1:length(t_nd)
       uVals(k,:) = (-Kr*dr_nd(k,:).' - P*dv_nd(k,:).').'; 
    end
    
    % Dimensionalize the relative quantities:
    dr_dim = dr_nd * L_scale;    % [km]
    dv_dim = dv_nd * V_scale;      % [km/s]
    u_dim  = uVals * A_scale;      % [km/s^2]
    
    % Apply saturation on the dimensional control input if flag is set:
    if flag_saturate
        for k = 1:length(t_nd)
            norm_u = norm(u_dim(k,:));
            if norm_u > u_sat_dim
                u_dim(k,:) = u_dim(k,:) * (u_sat_dim / norm_u);
            end
        end
    end
    
    % Compute magnitudes (dimensional):
    dr_mag_dim = sqrt(sum(dr_dim.^2,2));
    dv_mag_dim = sqrt(sum(dv_dim.^2,2));
    u_mag_dim  = sqrt(sum(u_dim.^2,2));
    
    % Convert nondimensional time to dimensional time in hours:
    t_hrs = (t_nd * T_scale) / 3600;
    
    % ---------- Relative State and Control Plots (Dimensional) ----------
    
    % Figure f1: Error States (delta r and delta v) vs time.
    f1 = figure('Position',[x_offset, y_offsets(1), figW, figH],'Color','w');
    subplot(2,1,1)
      plot(t_hrs, dr_dim(:,1), 'r','LineWidth',2); hold on;
      plot(t_hrs, dr_dim(:,2), 'g','LineWidth',2);
      plot(t_hrs, dr_dim(:,3), 'b','LineWidth',2);
      plot(t_hrs, dr_mag_dim, 'k--','LineWidth',2);
      xlabel('Time (hrs)','FontSize',20);
      ylabel('\delta r (km)','FontSize',20);
      title([scenarioTag,': \delta r Components & Magnitude'],'FontSize',20);
      legend('\delta r_1','\delta r_2','\delta r_3','|\delta r|','Location','best','FontSize',20);
      grid on; set(gca,'FontSize',20); hold off;
    subplot(2,1,2)
      plot(t_hrs, dv_dim(:,1), 'r','LineWidth',2); hold on;
      plot(t_hrs, dv_dim(:,2), 'g','LineWidth',2);
      plot(t_hrs, dv_dim(:,3), 'b','LineWidth',2);
      plot(t_hrs, dv_mag_dim, 'k--','LineWidth',2);
      xlabel('Time (hrs)','FontSize',20);
      ylabel('\delta v (km/s)','FontSize',20);
      title([scenarioTag,': \delta v Components & Magnitude'],'FontSize',20);
      legend('\delta v_1','\delta v_2','\delta v_3','|\delta v|','Location','best','FontSize',20);
      grid on; set(gca,'FontSize',20); hold off;
    
    % Figure f2: Control Input vs time.
    f2 = figure('Position',[x_offset, y_offsets(2), figW, figH],'Color','w');
      plot(t_hrs, u_dim(:,1), 'r','LineWidth',2); hold on;
      plot(t_hrs, u_dim(:,2), 'g','LineWidth',2);
      plot(t_hrs, u_dim(:,3), 'b','LineWidth',2);
      plot(t_hrs, u_mag_dim,      'k--','LineWidth',2);
      xlabel('Time (hrs)','FontSize',20);
      ylabel('Control Input (km/s^2)','FontSize',20);
      title([scenarioTag,': Control Input & Magnitude'],'FontSize',20);
      legend('u_1','u_2','u_3','|u|','Location','best','FontSize',20);
      grid on; set(gca,'FontSize',20); hold off;
    
    % Figure f3: 3D Plot of Relative Error State.
    f3 = figure('Position',[x_offset, y_offsets(3), figW, figH],'Color','w');
      plot3(dr_dim(:,1), dr_dim(:,2), dr_dim(:,3),'LineWidth',2);
      xlabel('\delta r_1 (km)','FontSize',20);
      ylabel('\delta r_2 (km)','FontSize',20);
      zlabel('\delta r_3 (km)','FontSize',20);
      title([scenarioTag,': 3D Plot of \delta r(t)'],'FontSize',20);
      grid on; axis equal; set(gca,'FontSize',20);
    
    % Figure f4: Plot Δ_i for i = 1,2,3 and Overall Δ.
    Delta1 = sqrt( (dr_dim(:,1)).^2 + (dv_dim(:,1)).^2 );
    Delta2 = sqrt( (dr_dim(:,2)).^2 + (dv_dim(:,2)).^2 );
    Delta3 = sqrt( (dr_dim(:,3)).^2 + (dv_dim(:,3)).^2 );
    Delta_total = sqrt( (dr_mag_dim).^2 + (dv_mag_dim).^2 );
    
    f4 = figure('Position',[x_offset, y_offsets(4), figW, figH],'Color','w');
      plot(t_hrs, Delta1, 'r', 'LineWidth',2); hold on;
      plot(t_hrs, Delta2, 'g', 'LineWidth',2);
      plot(t_hrs, Delta3, 'b', 'LineWidth',2);
      plot(t_hrs, Delta_total, 'k--', 'LineWidth',2);
      xlabel('Time (hrs)','FontSize',20);
      ylabel('\Delta (km, km/s)','FontSize',20);
      title([scenarioTag,': \Delta_i = sqrt(\delta r_i^2 + \delta v_i^2) and Overall \Delta'],'FontSize',20);
      legend('\Delta_1','\Delta_2','\Delta_3','Overall \Delta','Location','best','FontSize',20);
      grid on; set(gca,'FontSize',20); hold off;
    
    % Figure f5: Absolute Orbit Propagation and 3D Plot (High Resolution)
    t_nd_high = linspace(tspan_nd(1), tspan_nd(2), 1000);  % High-res time vector
    dr_nd_high = interp1(t_nd, dr_nd, t_nd_high);          % Interpolate relative state
    rChief_nd_high = propagateOrbit_nd(1, e_chief, i_chief, Omega_chief, omega_chief, nu_chief, 1, t_nd_high);
    rDeputy_nd_high = rChief_nd_high + dr_nd_high;
    rChief_dim = rChief_nd_high * L_scale;
    rDeputy_dim = rDeputy_nd_high * L_scale;
    
    f5 = figure('Position',[x_offset, y_offsets(5), figW, figH],'Color','w');
      plot3(rChief_dim(:,1), rChief_dim(:,2), rChief_dim(:,3),'b','LineWidth',2); hold on;
      plot3(rDeputy_dim(:,1), rDeputy_dim(:,2), rDeputy_dim(:,3),'r:','LineWidth',2);
      plot3(rChief_dim(end,1), rChief_dim(end,2), rChief_dim(end,3),'bo','MarkerSize',20,'MarkerFaceColor','b');
      plot3(rDeputy_dim(end,1), rDeputy_dim(end,2), rDeputy_dim(end,3),'ro','MarkerSize',10,'MarkerFaceColor','r');
      xlabel('X (km)','FontSize',20);
      ylabel('Y (km)','FontSize',20);
      zlabel('Z (km)','FontSize',20);
      title([scenarioTag,': 3D Chief and Deputy Orbits'],'FontSize',20);
      legend('Chief Orbit','Deputy Orbit','Chief Final','Deputy Final','Location','best','FontSize',20);
      grid on; axis equal; set(gca,'FontSize',20);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLOSED-LOOP DYNAMICS (nondimensional):
function dxdt = closedLoopDynamics(x, Kr, P)
    dr = x(1:3);
    dv = x(4:6);
    dxdt = [dv; -Kr * dr - P * dv];
end

% OE2RV: Convert orbital elements to Cartesian ECI state vectors (dimensional):
function [rECI, vECI] = oe2rv(a, e, i_deg, Omega_deg, omega_deg, nu_deg, mu)
    i = deg2rad(i_deg);
    Omega = deg2rad(Omega_deg);
    omega = deg2rad(omega_deg);
    nu = deg2rad(nu_deg);
    p = a * (1 - e^2);
    r_orb = p / (1 + e*cos(nu));
    r_pf = [r_orb*cos(nu); r_orb*sin(nu); 0];
    h = sqrt(mu * p);
    v_pf = [-mu/h*sin(nu); mu/h*(e+cos(nu)); 0];
    R = [ cos(Omega)*cos(omega)-sin(Omega)*sin(omega)*cos(i), -cos(Omega)*sin(omega)-sin(Omega)*cos(omega)*cos(i),  sin(Omega)*sin(i);
          sin(Omega)*cos(omega)+cos(Omega)*sin(omega)*cos(i), -sin(Omega)*sin(omega)+cos(Omega)*cos(omega)*cos(i), -cos(Omega)*sin(i);
          sin(omega)*sin(i),                                cos(omega)*sin(i),                                 cos(i)];
    rECI = R * r_pf;
    vECI = R * v_pf;
end

% propagateOrbit_nd: Propagate orbit in nondimensional units.
function rECI_nd = propagateOrbit_nd(a, e, i_deg, Omega_deg, omega_deg, nu0_deg, mu, tvec)
    i = deg2rad(i_deg);
    Omega = deg2rad(Omega_deg);
    omega = deg2rad(omega_deg);
    nu0 = deg2rad(nu0_deg);
    n = sqrt(mu / a^3);
    E0 = 2*atan( sqrt((1-e)/(1+e)) * tan(nu0/2) );
    M0 = E0 - e*sin(E0);
    N = length(tvec);
    rECI_nd = zeros(N,3);
    R = [ cos(Omega)*cos(omega)-sin(Omega)*sin(omega)*cos(i), -cos(Omega)*sin(omega)-sin(Omega)*cos(omega)*cos(i),  sin(Omega)*sin(i);
          sin(Omega)*cos(omega)+cos(Omega)*sin(omega)*cos(i), -sin(Omega)*sin(omega)+cos(Omega)*cos(omega)*cos(i), -cos(Omega)*sin(i);
          sin(omega)*sin(i),                                cos(omega)*sin(i),                                 cos(i)];
    for j = 1:N
        M = M0 + n*tvec(j);
        E = M;
        for k = 1:50
            E_new = M + e*sin(E);
            if abs(E_new - E) < 1e-8, break; end
            E = E_new;
        end
        nu = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
        r = a * (1 - e*cos(E));
        r_pf = [r*cos(nu); r*sin(nu); 0];
        rECI_nd(j,:) = (R * r_pf).';
    end
end












