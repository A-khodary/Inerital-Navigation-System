function [p_INS, v_INS] = INS(p_ant, v_ant,qhat, a, g,T)
% d2r =pi/180;
% 
% delta_ang = 16.000039000368773*d2r;
% eangles  = quat2euler(qhat);
% eangles.psi = eangles.psi - delta_ang;
% qhat = euler2quat(eangles);


lambda = p_ant(1);
h = p_ant(3);

vn = v_ant(1);
ve = v_ant(2);
vd = v_ant(3);

Rlambda = Rlambda_m(lambda);
Rphi = Rphi_m(lambda);

dv_INS = rot(qhat)'*a - [0; 0; g] + [ ((-ve^2*sin(lambda))/((Rphi + h)*cos(lambda)) + (vn*vd)/(Rlambda + h) );
                                      ((ve*vn*sin(lambda))/((Rphi + h)*cos(lambda)) + (ve*vd)/(Rphi + h) ); 
                                      ((-ve^2)/(Rphi + h) - (vn^2)/(Rlambda + h) )];

                                  
v_INS = dv_INS*T + v_ant;

%v_INS = (rot(qhat)'*a - [0; 0; g])*T + v_ant; % Cálculo das velocidade NED pelo acelerômetro

dp_INS = [1/(Rlambda + h) 0 0;
    0                 1/((Rphi + h)*cos(lambda)) 0;
    0                 0                       -1]*v_INS;
p_INS = dp_INS*T+p_ant;
end