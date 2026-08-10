function [D,G,I]=build_ops0(Grid) 
% author: your name 
% date: today 
% description: 
% This function computes the discrete divergence and gradient matrices on a 
% regular staggered grid using central difference approximations. The 
% discrete gradient assumes homogeneous boundary conditions. 
% Input: 
% Grid = structure containing all pertinent information about the grid. 
% Output: 
% D = Nx by Nx+1 discrete divergence matrix 
% G = Nx+1 by Nx discrete gradient matrix 
% I = Nx by Nx identity matrix 
% 
% Example call: 
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10; 
% >> Grid = build_grid(Grid); 
% >> [D,G,I]=build_ops(Grid);
Nx=Grid.Nx;

B=zeros(Nx+1,2);
B(:,1)=-1;
B(:,2)=1;
D=1/Grid.dx*spdiags(B,[0 1],Nx,Nx+1);

G=-D';
dof_f_bnd = [Grid.dof_f_xmin Grid.dof_f_xmax]; % all dof's on boundary
G(dof_f_bnd,:) = 0;

I=speye(Nx);

