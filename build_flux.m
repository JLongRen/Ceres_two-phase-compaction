function [Flux] = build_flux(Flux) 
Nt=Flux.Nt;
Nf=Flux.Nf;
Flux.h=zeros(Nt,Nf);
Flux.phi=zeros(Nt,Nf);
Flux.H=zeros(Nt,Nf);
Flux.Phi=zeros(Nt,Nf);

Flux.source=zeros(Nt,Nf-1);

Flux.wf=zeros(Nt,Nf);
Flux.wF=zeros(Nt,Nf);
end 



