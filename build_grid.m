function [Grid] = build_grid(Grid) 
% author: Jialong Ren 
% date: 10 2019 
% Description: 
% This function computes takes in minimal definition of the computational 
% domain and grid and computes all containing all pertinent information 
% about the grid. 
% 
% Input: 
% Grid.xmin = left boundary of the domain 
% Grid.xmax = right bondary of the domain 
% Grid.Nx = number of grid cells 
% 
% Output: (suggestions) 
% Grid.Lx = scalar length of the domain 
% Grid.dx = scalar cell width 
% Grid.Nfx = number of fluxes in x-direction 
% Grid.xc = Nx by 1 column vector of cell center locations 
% Grid.xf = Nfx by 1 column vector of cell face locations 
% Grid.dof = Nx by 1 column vector from 1 to N containig the degrees of freedom, i.e. cell numbers 
% Grid.dof_xmin = scalar cell degree of freedom corrsponding to the left boundary 
% Grid.dof_xmax = scalar cell degree of freedom corrsponding to the right boundary 
% Grid.dof_f_xmin = scalar face degree of freedom corrsponding to the left boundary 
% Grid.dof_f_xmax = scalar face degree of freedom corrsponding to the right boundary 
% 
% Example call: 
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10; 
% >> Grid = build_grid(Grid);
 Grid.Lx =Grid.xmax-Grid.xmin;

 Grid.dx = Grid.Lx/Grid.Nx;
if ~isfield(Grid,'psi_x0')
    Grid.psi_x0 = 'xmin_ymin'; 
end
if ~isfield(Grid,'psi_dir')
    Grid.psi_dir = 'xy'; 
end
 if ~isfield(Grid,'ymin')
     Grid.ymin = 0; 
     Grid.Ny=1;
 else
       Grid.Ly =Grid.ymax-Grid.ymin;
         Grid.dy = Grid.Ly/Grid.Ny;
 end
 if ~isfield(Grid,'ymax')
     Grid.ymax = Grid.dx; 
 end
 
 Grid.Ly =Grid.ymax-Grid.ymin;
 Grid.dy = Grid.Ly/Grid.Ny;
  
 Grid.N =Grid.Nx*Grid.Ny; 
 Grid.Nfx = (Grid.Nx+1)*Grid.Ny;
 Grid.Nfy = (Grid.Ny+1)*Grid.Nx;
 Grid.Nf=Grid.Nfx+Grid.Nfy;
 
 Xc = linspace(0.5*Grid.dx,(Grid.Nx-0.5)*Grid.dx,Grid.Nx);
 Yc = linspace(0.5*Grid.dy,(Grid.Ny-0.5)*Grid.dy,Grid.Ny);
%  [Grid.xc0,Grid.yc0] = meshgrid(Xc,Yc);
%  Grid.xc=reshape(Grid.xc0,Grid.N,1);
%  Grid.yc=reshape(Grid.yc0,Grid.N,1);
 Grid.xc =Grid.xmin+Xc;
 Grid.yc =Grid.ymin+Yc;
 
 Xf = linspace(0,Grid.Lx,Grid.Nx+1);
 %[Grid.xf1,Grid.yf1]=meshgrid(Xf,Yc);
%  Grid.xf(1:Grid.Nfx,1)=reshape(Grid.xf1,Grid.Nfx,1);  
%  Grid.yf(1:Grid.Nfx,1)=reshape(Grid.yf1,Grid.Nfx,1);
 Grid.xf=Grid.xmin+Xf;
 
 Yf = linspace(0,Grid.Ly,Grid.Ny+1);
%  [Grid.xf2,Grid.yf2]=meshgrid(Xc,Yf);
%  Grid.xf(Grid.Nfx+1:Grid.Nf,1)=reshape(Grid.xf2,Grid.Nfy,1);  
%  Grid.yf(Grid.Nfx+1:Grid.Nf,1)=reshape(Grid.yf2,Grid.Nfy,1);
  Grid.yf=Yf;

 Ntdof=(1:Grid.N);% temporary dof 1:Grid.N
 Grid.dof = reshape(Ntdof,Grid.Ny,Grid.Nx);
 
 Ntfdof=(1:Grid.Nf);% temporary  dof 1:Grid.Nf
 Grid.dof_fx =reshape(Ntfdof(1:Grid.Nfx),Grid.Ny,Grid.Nx+1);
 Grid.dof_fy =reshape(Ntfdof(Grid.Nfx+1:Grid.Nf),Grid.Ny+1,Grid.Nx);

  Grid.dof_xmin = Grid.dof(:,1);
  Grid.dof_xmax =Grid.dof(:,Grid.Nx);
  Grid.dof_ymin = Grid.dof(1,:)';
  Grid.dof_ymax =Grid.dof(Grid.Ny,:)';
  
  Grid.dof_f_xmin = Grid.dof_fx(:,1);
  Grid.dof_f_xmax = Grid.dof_fx(:,Grid.Nx+1);
  Grid.dof_f_ymin = Grid.dof_fy(1,:)';
  Grid.dof_f_ymax = Grid.dof_fy(Grid.Ny+1,:)';

  Grid.dz=Grid.dx;
  if~isfield(Grid,'geo')
  Grid.V=Grid.dx*Grid.dy*Grid.dz*ones(Grid.N,1);
  Grid.A(1:Grid.Nfx)=Grid.dy*Grid.dz*ones(Grid.Nfx,1);
  Grid.A(Grid.Nfx+1:Grid.Nf)=Grid.dx*Grid.dz*ones(Grid.Nfy,1);

  elseif Grid.geo==2  %spherical
        Grid.V=4/3*pi*((Grid.xc+0.5*Grid.dx).^3-(Grid.xc-0.5*Grid.dx).^3)';
  Grid.A(1:Grid.Nfx)=4*pi*Grid.xf.^2;
  
  else
      
 Grid.dof_f_bnd=[ Grid.dof_f_xmin; Grid.dof_f_xmax; Grid.dof_f_ymin; Grid.dof_f_ymax ];
      error('geo type not defined')
  end
  Grid.A=Grid.A';
end 



