function [Opropeties] = set_Opropeties(Nt,Nx) 
Opropeties.wa=zeros(Nt,Nx);
Opropeties.Kappac=zeros(Nt,Nx); % thermal conductivity K
Opropeties.KdH=zeros(Nt,Nx+1);   %  KddTdH
Opropeties.Kdphi=zeros(Nt,Nx+1);   %KddTdphi
Opropeties.Cp=ones(Nt,Nx);
Opropeties.pp=zeros(Nt,Nx+1);
Opropeties.phi=zeros(Nt,Nx);
Opropeties.dphi=zeros(Nt,Nx);

%wa,Kappac,Cpc,ppc,Kappac2,x
end 



