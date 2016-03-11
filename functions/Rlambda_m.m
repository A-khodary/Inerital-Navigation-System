function [Rlambda] = Rlambda_m(lambda);

e = 0.0818;        % excentricidade
a = 6378137;       % raio equatiorial m

Rlambda = (a*(1-e^2))/((1-e^2*sin(lambda)^2)^1.5);