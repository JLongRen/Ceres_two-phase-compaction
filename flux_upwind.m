function [A] = flux_upwind(q,Nx,Nfx) 
% author: Jialong REN
% date: 2019 Nov 26th
% Description: 
% This function computes the upwind flux matrix from the flux vector. 
% 
% Input: 
% q = Nf by 1 flux vector from the flow problem. 
% Grid = structure containing all pertinent information about the grid. 
% % Output: 
% A = Nf by Nf matrix contining the upwinded fluxes 
% Nx = Grid.Nx; 
% 
% Nfx = Grid.Nfx; % # of x faces Nfy = Grid.Nfy; % # of y faces

    %% One dimensinal 
    qx=q(1:Nfx);
        qn =min(qx(1:Nfx-1),0); 
        qp = max(qx(2:Nfx),0);  
        A = spdiags(qn,0,Nfx,Nx)+spdiags(qp,-1,Nfx,Nx);

end 
