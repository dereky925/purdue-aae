function dx = linearizedTwoPendulumTwoCart(t,x,A,B)

    u = 0;

    dx = A * x + B * u;

end