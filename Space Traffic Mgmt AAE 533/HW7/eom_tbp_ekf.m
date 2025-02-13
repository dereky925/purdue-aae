function [dqdt] = eom_tbp_ekf(~,q,mu,M,Q)

x  = q(1);
y  = q(2);
z  = q(3);
dx = q(4);
dy = q(5);
dz = q(6);
P = reshape(q(7:42),6,6);

r2 = x*x + y*y + z*z;
r = sqrt(r2);

s = x/r;
t = y/r;
u = z/r;

f1 = mu/r2;
f2 = f1/r;

% Compute derivatives of the potential function (the acceleration)
Ux = -f1*s;
Uy = -f1*t;
Uz = -f1*u;

Gxx = f2*(3.0*s*s - 1.0);
Gxy = f2*(3.0*s*t);
Gxz = f2*(3.0*s*u);
Gyy = f2*(3.0*t*t - 1.0);
Gyz = f2*(3.0*t*u);
Gzz = f2*(3.0*u*u - 1.0);

G = [Gxx, Gxy, Gxz; Gxy, Gyy, Gyz; Gxz, Gyz, Gzz];

% Assemble the dynamics Jacobian
F = [zeros(3,3), eye(3); G, zeros(3)];

% Time rates for the mean and covariance
dmdt = [dx; dy; dz; Ux; Uy; Uz];
dPdt = F*P + P*F' + M*Q*M';

% Assemble time derivates for output
dqdt = [dmdt;dPdt(:)];

end






























