function [Lim,Lex] = build_ade_ops(Nx,Nfx,D,q,I,dt) 
% file: build_ade_ops.m 
% author: Jialong Ren 
% date: 1206 2019 
% Description: 
% Assembles the function handles for the implicit and explicit matrices 
% arising from the discretization of the advection-dispersion equation. 
% The hydrodynamic dispersion term is always treated implicitly. The 
% advection term is treated using the theta-method: 
% theta = 1 is the explicit Forward Euler method, 
% theta = 0 is thicit Backward Euler method.

% Input: 
% Grid = grid data structure. 
% I = Nc by Nc identity matrix. 
% D = Nc by Nf divergence matrix. 
% G = Nf by Nc gradient matrix. 
% q = Nf by 1 flux vector 
% phi = Nc by 1 vector of porosities. 
% Disp = Hydrodynamic dispersion. 
% dt = time step 
% Output: 
% Lim = @(theta) Implicit N by N matrix from the discretization of the 
% advection-diffusion equation. Lim should be an anonymous function 1 of theta. 
% Lex = @(theta) Explicit N by N matrix from the discretization of the 
% advection-diffusion equation. Lex should be an anonymous function of theta.
 A = flux_upwind(q,Nx,Nfx);
Lim = @(theta) I+theta*D*A*dt;
Lex = @(theta) I-(1-theta)*D*A*dt;