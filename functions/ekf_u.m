function [x,P] = ekf_u(xp, Pp, z, sys,opt)

H = sys.H;
R = sys.R;

m = size(xp,1);
epsilon = 10^(-1); % To prevent inversion of a singular matrix
K = Pp*H'*inv(H*Pp*H' + R + epsilon*eye(size(R)));
P = (eye(m) - K*H)*Pp;
switch opt
    case 0 
        x = xp + K*(z-H*xp);
    case 1
        x = xp + K*z; % z = deltaz;
    otherwise
        disp('Unknow method ekf');
end