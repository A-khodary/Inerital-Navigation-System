function [sys] = sysEuler_f(rates, xhat, varQ, varR, g, tau, T)

 E = diag(varQ);
         

 R = diag(varR);

% E = [ 0.0001 0      0    0      0       0;
%     0        0.0001 0    0      0       0;
%     0        0      0.01 0      0       0;
%     0        0      0    0.0001 0       0;
%     0        0      0    0      0.0001  0;
%     0        0      0    0      0       0.0001];

G = eye(6);

Q = G*E*G'*T;


% R = [  6 0 0 ;
%        0 6 0 ; 
%        0 0 6 ];


H = [ 1 0 0 0 0 0;
      0 1 0 0 0 0;
      0 0 1 0 0 0];


A = zeros(3, 3);

phi   = xhat(1);
theta = xhat(2);

bp = xhat(4);
bq = xhat(5);
br = xhat(6);

ps = rates(1);
qs = rates(2);
rs = rates(3);

A(1,1) = qs*cos(phi)*tan(theta)   - rs*sin(phi)*tan(theta)   - bq*cos(phi)*tan(theta)   + br*sin(phi)*tan(theta);
A(1,2) = qs*sin(phi)*sec(theta)^2 + rs*cos(phi)*sec(theta)^2 - bq*sin(phi)*sec(theta)^2 - br*cos(phi)*sec(theta)^2;
A(1,3) = 0;
A(1,4) = -1;
A(1,5) = sin(phi)*tan(theta);
A(1,6) = cos(phi)*tan(theta);

A(2,1) = - qs*sin(phi) - rs*cos(phi) + bq*cos(phi) - br*cos(phi);
A(2,2) = 0;
A(2,3) = 0;
A(2,4) = 0;
A(2,5) = - cos(phi);
A(2,6) = - sin(phi);


A(3,1) = qs*cos(phi)*sec(theta) - rs*sin(phi)*sec(theta) - bq*cos(phi)*sec(theta) + br*sin(phi)*sec(theta);
A(3,2) = qs*sin(phi)*sec(theta)*tan(theta) + rs*cos(phi)*sec(theta)*tan(theta) - bq*sin(phi)*sec(theta)*tan(theta) - br*cos(phi)*sec(theta)*tan(theta);
A(3,3) = 0;
A(3,4) = 0;
A(3,5) = - sin(phi)*sec(theta);
A(3,6) = - cos(phi)*sec(theta);

A(4,1) = 0;
A(4,2) = 0;
A(4,3) = 0;
A(4,4) = -1/tau(1);
A(4,5) = 0;
A(4,6) = 0;


A(5,1) = 0;
A(5,2) = 0;
A(5,3) = 0;
A(5,4) = 0;
A(5,5) = -1/tau(2);
A(5,6) = 0;

A(6,1) = 0;
A(6,2) = 0;
A(6,3) = 0;
A(6,4) = 0;
A(6,5) = 0;
A(6,6) = -1/tau(3);

Phi = eye(6) + A*T;

sys.Phi = Phi;
sys.G = G;
sys.H = H;
sys.Q = Q;
sys.Qeta = E;
sys.R = R;

Ef = 0.001*ones(6,6);
Eg = 0.001*ones(6,6);
Eh = 0.001*ones(3,6);
Mf = eye(6);
Mh = eye(3);
sys.alfa = 2000; % Sayed
sys.Ef = Ef;
sys.Eg = Eg;
sys.Eh = Eh;
% sys.M = M;
sys.Mf = Mf;
sys.Mh = Mh;

sys.gamma = 100000;