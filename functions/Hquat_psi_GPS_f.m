function [H] = Hquat_psi_GPS_f(q, g,opt)


q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);

r = rot(q);

switch opt
    case 1
        H = [eye(4), zeros(4,3)]; 
        
    case 0
        H = [eye(4), zeros(4,3)];
        H(5,1) = (2*q3*r(1,1) - 2*q0*r(1,2))/(r(1,1)^2);
        H(5,2) = (2*q2*r(1,1) - 2*q1*r(1,2))/(r(1,1)^2);
        H(5,3) = (2*q1*r(1,1) + 2*q2*r(1,2))/(r(1,1)^2);
        H(5,4) = (2*q0*r(1,1) + 2*q3*r(1,2))/(r(1,1)^2);
        H(5,5:7) = zeros(1,3);
end