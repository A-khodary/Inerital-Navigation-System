function [phi, theta, psi, bg] = EulerEKF_bg(z, rates, dt)
%
%
persistent H Q R
persistent x P
persistent tauR
persistent firstRun


if isempty(firstRun)
  H = [ 1 0 0 0 0 0;
        0 1 0 0 0 0;
        0 0 1 0 0 0];
    
  Q = [ 0.0001 0      0   0 0 0;
        0      0.0001 0   0 0 0;
        0      0      0.01 0 0 0;
        0      0      0   0.0001 0 0;
        0      0      0   0      0.0001  0;
        0      0      0   0      0       0.0001];
     
%   R = [  6 0 0 0;
%          0 6 0 0; 
%          0 0 6 0;
%          0 0 0 0.06];
     
    R = [  6 0 0 ;
           0 6 0 ; 
           0 0 6 ];

  x = [0 0 0 0.01 0.01 0.01]';  
  P = 10*eye(6);
  
  
  firstRun = 1;  
end


A = Ajacob(x, rates, dt);

xp = fx(x, rates, dt);
Pp = A*P*A' + Q;

K = Pp*H'*inv(H*Pp*H' + R);

% x = xp + K*(z - H*xp);
x = xp + K*normalize_angle_f(z-H*xp,-pi);

P = Pp - K*H*Pp;


phi   = x(1);
theta = x(2);
psi   = x(3);

bg = x(4:6);


%------------------------------
function xp = fx(xhat, rates, dt)
%
%
tauR(1,1) = 626.8115;
tauR(2,1) = 6468.0515;
tauR(3,1) = 602.5784; 

phi   = xhat(1);
theta = xhat(2);

ps = rates(1);
qs = rates(2);
rs = rates(3);

bp = xhat(4);
bq = xhat(5);
br = xhat(6);


xdot = zeros(3, 1);
xdot(1) = ps + qs*sin(phi)*tan(theta) + rs*cos(phi)*tan(theta) -bp - bq*sin(phi)*tan(theta) - br*cos(phi)*tan(theta);
xdot(2) =     qs*cos(phi)            - rs*sin(phi) - bq*cos(phi) - br*sin(phi);
xdot(3) =     qs*sin(phi)*sec(theta) + rs*cos(phi)*sec(theta) - bq*sin(phi)*sec(theta) - br*cos(phi)*sec(theta);
xdot(4) = -1/tauR(1)*bp;
xdot(5) = -1/tauR(2)*bq;
xdot(6) = -1/tauR(3)*br;
xp = xhat + xdot*dt;


%------------------------------
function A = Ajacob(xhat, rates, dt)
%
%
A = zeros(3, 3);

tauR(1,1) = 626.8115;
tauR(2,1) = 6468.0515;
tauR(3,1) = 602.5784; 


phi   = xhat(1);
theta = xhat(2);

bp = xhat(4);
bq = xhat(5);
br = xhat(6);

ps = rates(1);
qs = rates(2);
rs = rates(3);

A(1,1) = qs*cos(phi)*tan(theta)   - rs*sin(phi)*tan(theta) -bq*cos(phi)*tan(theta) + br*sin(phi)*tan(theta);
A(1,2) = qs*sin(phi)*sec(theta)^2 + rs*cos(phi)*sec(theta)^2 -bq*sin(phi)*sec(theta)^2 -br*cos(phi)*sec(theta)^2;
A(1,3) = 0;
A(1,4) = -1;
A(1,5) = sin(phi)*tan(theta);
A(1,6) = cos(phi)*tan(theta);

A(2,1) = -qs*sin(phi) - rs*cos(phi) + bq*cos(phi)-br*cos(phi);
A(2,2) = 0;
A(2,3) = 0;
A(2,4) = 0;
A(2,5) = -cos(phi);
A(2,6) = -sin(phi);


A(3,1) = qs*cos(phi)*sec(theta) - rs*sin(phi)*sec(theta) - bq*cos(phi)*sec(theta) + br*sin(phi)*sec(theta);
A(3,2) = qs*sin(phi)*sec(theta)*tan(theta) + rs*cos(phi)*sec(theta)*tan(theta) - bq*sin(phi)*sec(theta)*tan(theta) - br*cos(phi)*sec(theta)*tan(theta);
A(3,3) = 0;
A(3,4) = 0;
A(3,5) = -sin(phi)*sec(theta);
A(3,6) = -cos(phi)*sec(theta);

A(4,1) = 0;
A(4,2) = 0;
A(4,3) = 0;
A(4,4) = -1/tauR(1);
A(4,5) = 0;
A(4,6) = 0;


A(5,1) = 0;
A(5,2) = 0;
A(5,3) = 0;
A(5,4) = 0;
A(5,5) = -1/tauR(2);
A(5,6) = 0;

A(6,1) = 0;
A(6,2) = 0;
A(6,3) = 0;
A(6,4) = 0;
A(6,5) = 0;
A(6,6) = -1/tauR(3);

A = eye(6) + A*dt;
