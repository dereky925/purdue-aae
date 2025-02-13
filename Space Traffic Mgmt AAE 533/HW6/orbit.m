
function f = orbit(t, x)

    mu = 3.986004418e5;
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    v1 = x(4);
    v2 = x(5);
    v3 = x(6);
    normx = sqrt(x1^2+x2^2+x3^2);
    
    vdot1 = ((x1)/(normx^3))*(-mu);
    vdot2 = ((x2)/(normx^3))*(-mu);
    vdot3 = ((x3)/(normx^3))*(-mu);
    xdot1 = v1;
    xdot2 = v2;
    xdot3 = v3;
    
    f = [xdot1; xdot2; xdot3; vdot1; vdot2; vdot3];
    
end
    
    