function [Omega] = Omegacalc(omega)
% [Omega] = Omegacalc(omega)
% Omega = [0           -omega(1)    -omega(2)    -omega(3);
%          omega(1)     0            omega(3)    -omega(2);
%          omega(2)    -omega(3)     0            omega(1);
%          omega(3)     omega(2)    -omega(1)     0];

Omega = [0           -omega(1)    -omega(2)    -omega(3);
         omega(1)     0            omega(3)    -omega(2);
         omega(2)    -omega(3)     0            omega(1);
         omega(3)     omega(2)    -omega(1)     0];