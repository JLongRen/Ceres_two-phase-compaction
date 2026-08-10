function [D,G,I]=build_ops(Grid) 
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
D2=1/Grid.dx*spdiags(B,[0 1],Nx,Nx+1);
G=-D2';
G(1,1)=0;
G(Nx+1,Nx)=0;
R=spdiags(Grid.xf.^2',0,Nx+1,Nx+1);
S=spdiags(1./(Grid.xc.^2'),0,Nx,Nx);
D=S*D2*R;
I=speye(Nx);

