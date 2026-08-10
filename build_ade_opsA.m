function [Lim,Lex] = build_ade_opsA(D,A,I,dt) 

Lim = @(theta) I+theta*D*A*dt;
Lex = @(theta) I-(1-theta)*D*A*dt;