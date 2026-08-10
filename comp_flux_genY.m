function [q] = comp_flux_genY(flux,res,u,Grid,BC,n1,n2)
% author: Marc Hesse
% date: 22 Feb 2019
% Description:
% Computes the fuxes on the interior from flux(u) and reconstructs the
% fluxes on the boundary faces from the residuals in the adjacent boundary
% cells, res(u). 

% Input:
% flux = anonymous function computing the flux (correct in the interior)
% res = anonymous function computing the residual 
% u = vector of 'flux potential' (head, temperature,electric field,...)
% Grid = structure containing pertinent information about the grid
% BC = structure containing pertinent information about BC's
% 
% Output:
% q = correct flux everywhere
%
% Example call:


%% Compute interior fluxes
q=flux(u);

%% Compute boundary fluxes
% 1) Identify the faces and cells on the boundary
dof_cell =[BC.dof_dir; BC.dof_neu];
dof_face = [BC.dof_f_dir; BC.dof_f_neu];
% 2) Determine sign of flux: Convention is that flux is positive in
%    coordinate direction. So the boundary flux, qb is not equal to q*n,
%    were n is the outward normal! 

 res_bnd =res(u,dof_cell-n1+1); 
sign = ones(size(dof_face,1),1)-2*ismember(dof_face,n2+1);

% 3) Compute residuals and convert them to bnd fluxes
q(dof_face-n1+1) = sign.*res_bnd.*Grid.V(dof_cell)./Grid.A(dof_face);
