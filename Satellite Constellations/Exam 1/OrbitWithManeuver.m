clc; clear; close all;

% Simulate an in-plane (in direction of velocity to be exact) maneuver
% to create an argument of perigee discontinuity

% -----------------------------
% 1. SETUP CONSTANTS & OPTIONS
% -----------------------------
mu  = MU;  % [km^3/s^2] Earth GM
Re  = EARTH_RADIUS;      % [km] Earth radius

% Maneuver time & magnitude (impulsive, in-plane)
tMan = 6000;          % [s] 10 min
dV   = 10.00;          % [km/s]

% -----------------------------
% 2. INITIAL ORBIT DEFINITION
% -----------------------------
% [a, e, i, RAAN, argPer, M0], angles in radians
% Example: a=7000 km, e=0.1, i=50 deg, RAAN=10 deg, argPer=90 deg, M0=0 deg

r_ijk = [656.0758  6596.8028 1473.0872]'; 
v_ijk = [-4.9631 -0.8127 5.7866]';
Y0 = [r_ijk,v_ijk];

[a, e, incl, O, w, nu] = eci2keplerian(r_ijk,v_ijk);

[r,v] = keplerian2eci(a+500, 0.0015, 15, 0, 100, nu);

Y0 = [r,v];

% -------------------------
% 3. PROPAGATE BEFORE MANEUVER
% -------------------------
tspan1  = [0, tMan];
options = odeset('RelTol',1e-12, 'AbsTol',1e-12);
[T1, Y1] = ode45(@(t, y) twoBodyJ2Dynamics(t, y, mu, J2, Re), tspan1, Y0, options);

% State at maneuver time
Yman = Y1(end, :)';

% ------------------------------------------------------
% 4. APPLY SMALL IN-PLANE DELTA-V (TANGENTIAL DIRECTION)
% ------------------------------------------------------
rVec = Yman(1:3);
vVec = Yman(4:6);
vHat = vVec / norm(vVec);   % unit velocity direction in-plane
Yman(4:6) = vVec + dV * vHat;  % apply impulse

% ------------------------
% 5. PROPAGATE AFTER MANEUVER
% ------------------------
tspan2 = [tMan, 6600];  % up to ~100 min total
[T2, Y2] = ode45(@(t, y) twoBodyJ2Dynamics(t, y, mu, J2, Re), tspan2, Yman, options);

% Combine time & states into single arrays
T = [T1; T2(2:end)];
Y = [Y1; Y2(2:end, :)];

% -------------------------
% 6. 3D ORBIT PLOT
% -------------------------
orbFig = OrbitPlotSetup();
plot3(Y(:,1), Y(:,2), Y(:,3),'b-','LineWidth',1.5);
orbFig.Position = [0 0 1000 1000];
% Mark the maneuver point
plot3(Y1(end,1), Y1(end,2), Y1(end,3), 'ro','MarkerSize',8,'LineWidth',2);

% -------------------------------
% 7. COMPUTE ARGUMENT OF PERIGEE VS TIME
% -------------------------------
nPoints = length(T);
wDeg    = zeros(nPoints,1);

for i = 1:nPoints

    [~,~,~,rval,wVal,~] = eci2keplerian([Y(i,1) Y(i,2) Y(i,3)],[Y(i,4) Y(i,5) Y(i,6)]);

    if rval > 0
        rval = rval-360;
    end
    RAAN(i) = rval;
    wDeg(i) = wVal;
end

% -------------------------
% 8. PLOT ARGUMENT OF PERIGEE VS. TIME
% -------------------------
figure('Name','Argument of Perigee','Color','w','Position',[1400 0 1400 1000]);
subplot(2,1,1)
hold on;grid minor
plot(T/60, wDeg, 'b-','LineWidth',2);
xlabel('Time [min]');
ylabel('\omega [deg]');
title('Argument of Perigee vs. Time (J_2 + Maneuver)');

% Indicate maneuver time
yLimits = ylim;
plot([tMan/60, tMan/60], [yLimits(1), yLimits(2)], 'r--','LineWidth',2);
legend('\omega(t)', 'Maneuver Time','Location','best');

subplot(2,1,2)
hold on;grid minor
plot(T/60, RAAN, 'b-','LineWidth',2);
xlabel('Time [min]');
ylabel('RAAN [deg]');
title('RAAN vs. Time (J_2 + Maneuver)');

% Indicate maneuver time
yLimits = ylim;
plot([tMan/60, tMan/60], [yLimits(1), yLimits(2)], 'r--','LineWidth',2);
legend('\omega(t)', 'Maneuver Time','Location','best');

subplotFont


% =========================================================================
% TWO-BODY + J2 DYNAMICS
% =========================================================================
function dydt = twoBodyJ2Dynamics(~, y, mu, J2, Re)
    % State y: [rx ry rz vx vy vz]
    rVec = y(1:3);
    vVec = y(4:6);
    r    = norm(rVec);

    % Two-body (point-mass) acceleration
    a2b  = -mu * rVec / r^3;

    % J2 perturbation
    z2   = rVec(3)^2;
    r2   = r^2;
    fac  = 1.5 * J2 * (mu * Re^2) / (r^5);
    aJ2x = fac * rVec(1) * (5*z2/r2 - 1);
    aJ2y = fac * rVec(2) * (5*z2/r2 - 1);
    aJ2z = fac * rVec(3) * (5*z2/r2 - 3);

    aJ2  = [aJ2x; aJ2y; aJ2z];
    aTot = a2b + aJ2;

    dydt = [vVec; aTot];
end

% =========================================================================
% ROTATION MATRIX about x=1, y=2, z=3 by angle th
% =========================================================================
function R = rotMatrix(axis, th)
    c = cos(th); s = sin(th);
    switch axis
        case 1
            R = [1  0  0
                 0  c -s
                 0  s  c];
        case 2
            R = [ c  0  s
                  0  1  0
                 -s  0  c];
        case 3
            R = [ c -s  0
                  s  c  0
                  0  0  1];
        otherwise
            R = eye(3);
    end
end

% =========================================================================
% SIMPLE KEPLER SOLVER: M = E - e*sin(E)
% =========================================================================
function E = solveKepler(M, e)
    E = M; 
    for k = 1:100
        f  = E - e*sin(E) - M;
        fp = 1 - e*cos(E);
        dE = -f/fp;
        E  = E + dE;
        if abs(dE) < 1e-14
            break;
        end
    end
end