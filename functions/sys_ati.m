function [sys] = sys_ati(sys,omega, q, sigma, R, g,m, tau, T, opt)
% q0 = q(1);
% q1 = q(2);
% q2 = q(3);
% q3 = q(4);
% 
% r = rot(q);

Omega = Omegacalc(omega);
Psi = Psicalc(q);

F = [0.5*Omega, -0.5*Psi;
    zeros(3,4)          , -diag(1./tau)];

% F = [0.5*Omega, -0.5*Psi, zeros(4,3);
%     zeros(3,4)          , -diag(1./tau(1:3)), zeros(3,3);
%     zeros(3,4)          , -diag(1./tau(4:6)), zeros(3,3)];

Phi = eye(7) + F*T;

%deltaPhi = 
omegap = sqrt(omega(1)^2 + omega(2)^2 + omega(3)^2);
Phi(1:4,1:4) = (eye(4)*cos(1/2*omegap*T) + Omega*inv(omegap)*sin(1/2*omegap*T));


G = [-0.5*Psi,  zeros(4,3);
          zeros(3,3),    eye(3,3)];

% G = [-0.5*Psi,  zeros(4,3), zeros(4,3);
%           zeros(3,3),    eye(3,3), zeros(3,3);
%           zeros(3,3),    zeros(3,3), eye(3,3)];

Gd = G*sqrt(T);
%Gd = G*T;

%Gd = (T*eye(10) + T^2*Phi)*G;
% H = [ -2*g*q2    2*g*q3     -2*g*q0      2*g*q1     0   0   0;
%        2*g*q1    2*g*q0      2*g*q3      2*g*q2     0   0   0;
%        g*q0   -g*q1     -g*q2      g*q3     0   0   0];%2*g*q0   -2*g*q1     -2*g*q2      2*g*q3, 4*g*q0      0            0        4*g*q3
% 
% H(4,1) = (2*q3*r(1,1) - 2*q0*r(1,2))/(r(1,1)^2 + r(1,2)^2);
% H(4,2) = (2*q2*r(1,1) - 2*q1*r(1,2))/(r(1,1)^2 + r(1,2)^2);
% H(4,3) = (2*q1*r(1,1) + 2*q2*r(1,2))/(r(1,1)^2 + r(1,2)^2);
% H(4,4) = (2*q0*r(1,1) + 2*q3*r(1,2))/(r(1,1)^2 + r(1,2)^2);
% H(4,5:7) = zeros(1,3);      
H = Hquat_m(q,g,m,opt);

Eeta1 = [sigma(1)^2, 0         , 0;
         0         , sigma(2)^2, 0;
         0         , 0         , sigma(3)^2];
    
Eeta2 = [sigma(4)^2, 0         , 0;
         0         , sigma(5)^2, 0;
         0         , 0         , sigma(6)^2];
     
% Eeta3 = [sigma(7)^2, 0         , 0;
%          0         , sigma(8)^2, 0;
%          0         , 0         , sigma(9)^2];
     
Eeta1eta2 = [Eeta1   , zeros(3);
             zeros(3), Eeta2];
% Eeta1eta2 = [Eeta1   , zeros(3), zeros(3);
%              zeros(3), Eeta2, zeros(3),
%              zeros(3), zeros(3), Eeta3];

Q = G*Eeta1eta2*G'*T;         
%Q = G*Eeta1eta2*G'*T^2;         
% if ~exist('sys.Q')
% Q = G*Eeta1eta2*G'*T;
% else
% Q = sys.Q;
% Q = Phi*Q*Phi' + G*Eeta1eta2*G';
% end

if opt  == 0
    R = R(1:3,1:3);
end
[p,n] = size(H);
[n,m] = size(G);

sys.Phi = Phi;
sys.G = G;
sys.H = H;
sys.Q = Q;
sys.Qeta = Eeta1eta2;
sys.R = R;
sys.Phi  = Phi;
sys.G    = G;
sys.Gd    = Gd;
sys.G_v  = zeros(n,p);
sys.H    = H;
sys.E    = eye(n);
sys.K_w  = zeros(p,m);
sys.K_v  = eye(p);
sys.J    = zeros(p,n);
sys.Q    = Q;
sys.S    = zeros(m,p);
% sys.Nf = 0.01*[0.5505    0.5505    0.5505    0.5505    0.1076    0.1076    0.1076 0 0 0]; sys.Ng_w = 0.01*[0.6721    0.6721    0.6721 0 0 0 0 0 0];
% sys.Nh = [0.0040    0.0281    0.0283    0.0288 0 0 0 0 0 0];
% sys.Nh = [0.0040    0.0284    0.0284    0.0284 0 0 0 0 0 0];
sys.Nf = (Phi - eye(7));
sys.Nf = [mean(abs(sys.Nf(1:4,1:7)))];
% sys.Nf = 100*[sys.Nf(1:4) 0.001*sys.Nf(5:7)];
% sys.Ng_w = 100*[mean(abs(Gd(1:4,1:3))) zeros(1,3)];
% sys.Nh = mean(H);
sys.Nf = 100*[sys.Nf(1:4) 0.001*sys.Nf(5:7)];
sys.Ng_w = 100*[mean(abs(Gd(1:4,1:3))) zeros(1,3)];
sys.Nh = mean(H);
%sys.Nh = [0.0040    0.0284    0.0284    0.0284 0 0 0 0 0 0];

sys.Nj   = zeros(1,n);
sys.Ng_v = zeros(1,p);
sys.Ne   = zeros(1,n);% 1e-12*ones(1,n);
sys.Nk_w = zeros(1,m);
sys.Nk_v = zeros(1,p);
sys.M1   = [1;1;1;1;0;0;0];
switch opt
    case 1
       %sys.M2   = [1;0;0;1;1;1];
       sys.M2   = [1;1;1;1;1;1];
    case 0
        sys.M2   = [1;1;1];
end

  sys.M1 = 0.001*sys.M1;
  sys.M2 = 0.001*sys.M2;

% sys.alfa = 5;
%sys.alfa = 4;
%sys.S_HInf = R;
sys.S_HInf = eye(7);

sys.L = eye(7);
% switch opt
%     case 1
%         sys.S_HInf = eye(6);
%     case 0
%         sys.S_HInf = eye(3);
% end
sys.theta = 0.1;
sys.gamma = 500000;
