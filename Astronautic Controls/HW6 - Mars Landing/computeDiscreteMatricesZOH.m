function [Acell, Bcell, ccell, xnext] = computeDiscreteMatricesZOH(fEoM, ...
                                     tspan, x0, uSeq, options)
%COMPUTEDISCRETEMATRICESZOH 
%  Computes (A_k, B_k, c_k) for each discrete step under ZOH
%  using numeric integration of an augmented system.
%
% Inputs:
%   fEoM   = handle to continuous dynamics: fEoM(t, x, u) -> xdot
%   tspan  = [t0, t1, ..., tN], time nodes
%   x0     = initial state at t0
%   uSeq   = cell array or matrix of control inputs [u0, u1, ..., u_{N-1}]
%   options= ODE solver options (optional)
%
% Outputs:
%   Acell{k} = A_k
%   Bcell{k} = B_k
%   ccell{k} = c_k
%   xnext{k} = x_{k+1} from the actual integration
%
% Each step integrates from t_k to t_{k+1} with u(t)=u_k (ZOH).

N = length(tspan) - 1;  % number of intervals
nx = length(x0);        % dimension of the state
nu = length(uSeq{1});   % dimension of the control

Acell = cell(1, N);
Bcell = cell(1, N);
ccell = cell(1, N);
xnext = cell(1, N);

xk = x0;  % current state
for k = 1:N
    
    tk     = tspan(k);
    tkp1   = tspan(k+1);
    uk     = uSeq{k};   % constant over [tk, tk+1]
    dt     = tkp1 - tk;
    
    % Build augmented initial condition:
    % Y = [x; reshape(Ix, nx^2,1); reshape(Iu, nx*nu,1)] 
    % where Ix=eye(nx), Iu= zeros(...) if you want partial w.r.t u, etc.
    
    Ix = eye(nx);
    Iu = zeros(nx, nu); 
    % or if your system is affine in u, you might keep partial derivatives in another structure
    Y0 = [xk; reshape(Ix, nx^2, 1); reshape(Iu, nx*nu, 1)];
    
    % Integrate from tk to tk+1:
    sol = ode45(@(t, Y) augEoM(t, Y, uk, fEoM, nx, nu), [tk tkp1], Y0, options);
    
    Yfinal = deval(sol, tkp1);  % Evaluate at end
    xkp1   = Yfinal(1:nx);      % This is x_{k+1}
    
    % Parse out partial derivatives wrt x:
    offset1 = nx; 
    PHI_x   = reshape(Yfinal(offset1+1 : offset1+nx^2), nx, nx);
    
    % partial derivatives wrt u (if you included them):
    offset2 = offset1 + nx^2;
    PHI_u   = reshape(Yfinal(offset2+1 : offset2+nx*nu), nx, nu);
    
    % Store results:
    xnext{k} = xkp1;
    
    % Now define A_k, B_k, c_k:
    % Suppose the continuous system is: dot{x} = A(t)*x + B(t)*u + c(t) (affine).
    % Then from the augmented integration, PHI_x ~ A_k, PHI_u ~ B_k, 
    % but there's also a constant shift (the integral of c(t)).
    
    Acell{k} = PHI_x;    % partial wrt x
    Bcell{k} = PHI_u;    % partial wrt u
    ccell{k} = xkp1 - PHI_x*xk - PHI_u*uk;
    
    % Update xk for next step
    xk = xkp1;
end

end

%% ---------------------------------------------------------------------
function dYdt = augEoM(t, Y, uk, fEoM, nx, nu)
% augEoM: builds an augmented system to track partial derivatives wrt x,u.
%
% Y = [x; vec(PHI_x); vec(PHI_u)]
% We'll differentiate each piece w.r.t time.

% Extract x from Y:
x    = Y(1:nx);
ix1  = nx;
ix2  = nx + nx^2;
ix3  = nx + nx^2 + nx*nu; 

% The continuous dynamics:
xdot = fEoM(t, x, uk);  % must be a function handle

% If fEoM is affine: xdot = A(t)*x + B(t)*u + c(t), 
% then partial_x = A(t), partial_u = B(t).
% Otherwise, you'd linearize or compute Jacobians numerically here.

% For demonstration, let's just do a generic approach:
A_t  = approximateJacobianX(@(xx) fEoM(t, xx, uk), x);   % Nx-by-Nx
B_t  = approximateJacobianU(@(uu) fEoM(t, x, uu), uk);   % Nx-by-Nu
% c(t) is included in xdot if the system is not purely linear

% Extract the PHI_x, PHI_u from Y:
PHI_x = reshape(Y(ix1+1:ix2), nx, nx);
PHI_u = reshape(Y(ix2+1:ix3), nx, nu);

% Now define derivatives:
dxdT    = xdot;
dPHI_xdT= A_t * PHI_x;  % chain rule
dPHI_udT= A_t * PHI_u + B_t;  % chain rule

% Pack up:
dYdt = [dxdT;
        reshape(dPHI_xdT, nx*nx, 1);
        reshape(dPHI_udT, nx*nu, 1)];
end

%% -------------- (example numeric Jacobian approximations) --------------
function A = approximateJacobianX(fun, x)
% Simple finite-difference for partial w.r.t x
nx = length(x);
A  = zeros(nx, nx);
epsVal = 1e-6;
f0 = fun(x);
for i = 1:nx
   dx      = zeros(nx,1);
   dx(i)   = epsVal;
   fPlus   = fun(x+dx);
   A(:, i) = (fPlus - f0)/epsVal;
end
end

function B = approximateJacobianU(fun, u)
% Simple finite-difference for partial w.r.t u
nu = length(u);
f0 = fun(u);
nx = length(f0);
B  = zeros(nx, nu);
epsVal = 1e-6;
for j = 1:nu
   du      = zeros(nu,1);
   du(j)   = epsVal;
   fPlus   = fun(u+du);
   B(:, j) = (fPlus - f0)/epsVal;
end
end