function [dxdt] = eom_robot(~, x, F)

dxdt = F*x;

end