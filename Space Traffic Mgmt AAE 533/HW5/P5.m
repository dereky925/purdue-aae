clear;clc;close all

x0 = [0 0 0 0.05 0.03 0.1]'; % [m, m/s]

% Seed
rng(1)

% Time
t0 = 0;
t = [t0:1:60]';

% Dynamics of Robot
F = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];

options = odeset('RelTol', 1e-10,'AbsTol',1e-15);
[~, X] = ode45(@eom_robot, t, x0, options, F);


% Create measurements
% STDev 1
z1 = X(1:3:30,1:3) + 1*randn(10, 3);
% STDev 3
z2 = X(31:3:60,1:3) + 3*randn(10, 3); 

z = [z1;z2];

% Initialize measurement vector and model matrix
Z = [];
H = [];
R = [];

% Iterate over each observation time
for k = 1:length(z)

    % Get time at kth observation
    tk = t(k)+1; % Not supposed to use full time vector...? it doesnt work
    
    % Concatenate measurements (assuming Z is a 3D matrix)
    Z = [Z; z(k,1); z(k,2); z(k,3)];
    
    % Concatenate model matrix
    Htilde = [1, 0, 0, 0, 0, 0; 
              0, 1, 0, 0, 0, 0;
              0, 0, 1, 0, 0, 0];
          
    Phik0 = [1, 0, 0, tk - t0, 0, 0; 
             0, 1, 0, 0, tk - t0, 0;
             0, 0, 1, 0, 0, tk - t0; 
             0, 0, 0, 1, 0, 0;
             0, 0, 0, 0, 1, 0;
             0, 0, 0, 0, 0, 1];
    
    % Concatenate the model matrix
    H = [H; Htilde * Phik0];

    if k < 11
        R = blkdiag(R, 1*eye(3)); 
    else
        R = blkdiag(R, 9*eye(3)); 
    end
end

% Least-squares estimate
xhat0 = (H' * H) \ (H' * Z);

L1 = xhat0(1) + xhat0(4)*[1:20];
L2 = xhat0(2) + xhat0(5)*[1:20];
L3 = xhat0(5) + xhat0(6)*[1:20];


% Matlab line fit validation
p = polyfit([1:20], z(:,1), 1);
slope = p(1);
intercept = p(2);
P = intercept + slope*[1:20];

% LUMVE
x_hat = inv([H'*inv(R)*H])*H'*inv(R)*Z
P_cov = inv([H'*inv(R)*H])

LUMVE1 = x_hat(1) + x_hat(4)*[1:20];
LUMVE2 = x_hat(2) + x_hat(5)*[1:20];
LUMVE3 = x_hat(3) + x_hat(6)*[1:20];

% Truth
truth1 = x0(4) * [1:3:60];
truth2 = x0(5) * [1:3:60];
truth3 = x0(6) * [1:3:60];

% Plotting =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
t_vec = 1:3:60;

a = figure;
subplot(3,1,1)
hold on
grid minor
scatter(NaN,NaN,70,'b','filled')
plot(NaN,NaN,'b--','LineWidth',1)
plot(NaN,NaN,'b:','LineWidth',4)
plot(NaN,NaN,'b','LineWidth',4)

scatter(t_vec,z(:,1),70,'b','filled') % Plot measurements
plot(t_vec,z(:,1),'b','LineWidth',1)
plot(t_vec,L1,'b--')                      % Plot Least squares
plot(t_vec,LUMVE1,'b:','LineWidth',2)     % Plot LUMVE
plot(t_vec,truth1,'b','LineWidth',2)      % Plot Truth

legend('X Measurements','X Least Squares','X LUMVE','X Truth','FontSize',20,'Location','best')
xlabel('Time [sec]','FontSize',20)
ylabel('Position [m]','FontSize',20)

% plot(P,'b--','LineWidth',3)

subplot(3,1,2)
hold on
grid minor
scatter(NaN,NaN,70,'r','filled')
plot(NaN,NaN,'r--','LineWidth',1)
plot(NaN,NaN,'r:','LineWidth',4)
plot(NaN,NaN,'r','LineWidth',4)
scatter(t_vec,z(:,2),70,'r','filled') % Plot measurements
plot(t_vec,z(:,2),'r','LineWidth',1)
plot(t_vec,L2,'r--')                      % Plot Least squares
plot(t_vec,LUMVE2,'r:','LineWidth',2)     % Plot LUMVE
plot(t_vec,truth2,'r','LineWidth',2)      % Plot Truth

legend('Y Measurements','Y Least Squares','Y LUMVE','Y Truth','FontSize',20,'Location','best')
xlabel('Time [sec]','FontSize',20)
ylabel('Position [m]','FontSize',20)

subplot(3,1,3)
hold on
grid minor
scatter(NaN,NaN,70,'g','filled')
plot(NaN,NaN,'g--','LineWidth',1)
plot(NaN,NaN,'g:','LineWidth',4)
plot(NaN,NaN,'g','LineWidth',4)
scatter(t_vec,z(:,3),70,'g','filled')  % Plot measurements
plot(t_vec,z(:,3),'g','LineWidth',1)
plot(t_vec,L3,'g--')                      % Plot Least squares
plot(t_vec,LUMVE3,'g:','LineWidth',2)     % Plot LUMVE
plot(t_vec,truth3,'g','LineWidth',2)      % Plot Truth

legend('Z Measurements','Z Least Squares','Z LUMVE','Z Truth','FontSize',20,'Location','best')
xlabel('Time [sec]','FontSize',20)
ylabel('Position [m]','FontSize',20)

a.Position = [100 100 1000 1800];








