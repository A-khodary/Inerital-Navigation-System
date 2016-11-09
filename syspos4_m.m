function [sys] = syspos4_m(p, v, q, varQ, varR, tau, T)

lambda  = p(1);
phi     = p(2);
h       = p(3);

vN = v(1);
vE = v(2);
vD = v(3);

Rlambda = Rlambda_m(lambda);
Rphi = Rphi_m(lambda);
RT = rot(q)';

F(1,1) = 0;
F(1,2) = 0;
F(1,3) = -(vN)/(Rlambda + h)^2;
F(2,1) =  (vE*sin(lambda))/((Rphi + h)*cos(lambda)^2);
F(2,2) = 0;
F(2,3) = -(vE)/((Rphi + h)^2*cos(lambda));
F(3,1) = 0; 
F(3,2) = 0;
F(3,3) = 0;

F(1,4) = inv(Rlambda + h);
F(1,5) = 0;
F(1,6) = 0;
F(2,4) = 0;
F(2,5) = inv((Rphi + h)*cos(lambda));
F(2,6) = 0;
F(3,4) = 0;
F(3,5) = 0;
F(3,6) = -1;

F(1:3, 7:9) = zeros(3,3);

F(4,1) = -(vE)^2/(Rphi + h) - (vE^2*sin(lambda)^2)/((Rphi + h)*cos(lambda)^2);
F(4,2) = 0;
F(4,3) = (vE^2*sin(lambda))/((Rphi + h)^2*cos(lambda)) -(vN*vD)/((Rlambda + h)^2);
F(5,1) = (vE*vN) /(Rphi + h) + (vE*vN*(sin(lambda)^2))/((Rphi + h)*(cos(lambda)^2));
F(5,2) = 0;
F(5,3) = -((vE*vN)*sin(lambda))/((Rphi + h)^2*cos(lambda)) - (vE*vD)/((Rphi + h)^2);
F(6,1) = 0;
F(6,2) = 0;
F(6,3) = (vE^2)/((Rphi + h)^2) + (vN^2)/((Rlambda + h)^2);

F(4,4) = vD/(Rlambda + h);
F(4,5) = -(2*vE*sin(lambda))/((Rphi + h)*cos(lambda));
F(4,6) = vN/(Rlambda + h);
F(5,4) = (vE*sin(lambda))/((Rphi + h)*cos(lambda));
F(5,5) = (vN*sin(lambda))/((Rphi + h)*cos(lambda)) + vD/(Rphi + h);
F(5,6) = vE/(Rphi + h);
F(6,4) = -(2*vN)/(h);
F(6,5) = -(2*vE)/(Rphi + h);
F(6,6) = 0;

F(4:6, 7:9) = -RT;

F(7:9, 1:6) = zeros(3,6);

F(7:9, 7:9) = -diag(1./tau);

Phi = 1*eye(9) + F*T;



L = [ zeros(6,3); RT];

G = [ zeros(3,6); 
      -RT, zeros(3,3); 
      zeros(3,3), eye(3,3)];

H = [ eye(3,3), zeros(3,3), zeros(3,3)];

        
E = diag(varQ);
         
Q = G*E*G'*T;

R = diag(varR);

Mf = 10^-3.* [ones(6,1); zeros(3,1)];
Mh = 10^-3 .* [ones(3,1)];

Nf = 10^2*sum(F(1:6,1:9) .* T) / 6;
Ng = 10^2*[sum(G(1:6,1:3))/6 , zeros(1,3)];
Nh = 0*sum(H(:,1:9))/3;

sys.Phi = Phi;
sys.L = L;
sys.G = G;
sys.H = H;
sys.Q = Q;
sys.Qeta = E;
sys.R = R;
sys.alfa = 1000; % Sayed
% sys.alfa = 0.1;
% sys.alfa = 1;
sys.Ef = Nf;
sys.Eg = Ng;
sys.Eh = Nh;
%sys.M = M;
sys.Mf = Mf;
sys.Mh = Mh;
             