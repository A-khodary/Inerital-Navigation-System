function [xp, Pp] = ekf_p (x,P,sys)
Q = sys.Q;
Phi = sys.Phi;
xp = Phi*x;
Pp = Phi*P*Phi' + Q;