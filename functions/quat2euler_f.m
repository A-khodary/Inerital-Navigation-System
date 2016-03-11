function [eangles] = quat2euler_f(q)

q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);


eangles.psi   = atan2((2*(q1*q2 + q0*q3)),(q0^2 + q1^2 - q2^2 - q3^2));
% sintheta = -2*(q1*q3 - q0*q2);
% eangles.theta = atan2(sintheta,2*(q1*q2 + q0*q3)*(1/sin(eangles.psi)));
eangles.theta = asin(-2*(q1*q3 - q0*q2));
% if norm(-2*(q1*q3 - q0*q2))>1
%     jj=0;
% end
eangles.phi   = atan2((2*(q2*q3 + q0*q1)),(q0^2 - q1^2 - q2^2 + q3^2));