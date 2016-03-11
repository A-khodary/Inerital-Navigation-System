function [p_IMU, v_IMU] = GPS_IMU(p_ant, v_ant,qhat, a, g,T)

% d2r =pi/180;
% 
% delta_ang =  16.000039000368773*d2r;
% eangles  = quat2euler(qhat);
% eangles.psi = eangles.psi - delta_ang;
% qhat = euler2quat(eangles);

lambda = p_ant(1);
h = p_ant(3);
v_IMU = (rot(qhat)'*a - [0; 0; g])*T + v_ant; % Cálculo das velocidade NED pelo acelerômetro

Rlambda = Rlambda_m(lambda);
Rphi = Rphi_m(lambda);

dp_IMU = [1/(Rlambda + h) 0 0;
    0                 1/((Rphi + h)*cos(lambda)) 0;
    0                 0                       -1]*v_IMU;
p_IMU = dp_IMU*T+p_ant;
