function [uc,T,Opropeties,Flux,ti]=wout5(R,To,t1,Nx,phi,th0,w,Dt,paraphysic,phis)
% R: radius 
%To: initail tmeperature , is a vector
%phi is weight fraction
rho_r=paraphysic.rho_r;
rho_i=paraphysic.rho_i;
rho_w=paraphysic.rho_w;
Cpr=paraphysic.Cpr;
Cpi=paraphysic.Cpi;
Cpw=paraphysic.Cpw;
k_i=paraphysic.k_i;
k_r=paraphysic.k_sp;
k_w=paraphysic.k_w;
Lh=paraphysic.Lh;

Grid.xmin = 0; Grid.xmax =R; Grid.Nx =Nx; %
Grid.geo=2;
Grid = build_grid(Grid); 
[D,G,I]=build_ops(Grid); 
Tas=273.15;

 Ta=th0; %ambient T
 yr2s=365.25*24*3600;
  Lt0=size(t1,2);
 Lt=floor(size(t1,2)/Dt);
 Flux.Nt=Lt;
  Flux.Nf=Nx+1;
 Flux=build_flux(Flux) ;
  
 Opropeties=set_Opropeties(Lt,Nx) ;
 wa=Opropeties.wa;
 Kappac=Opropeties.Kappac;  %thermal conductivity on face
Cpc=Opropeties.Cp;   %Capacity
ppc=Opropeties.pp;  %mean index

 phic=zeros(size(wa));
 Param.dof_dir = Grid.dof_xmax; % identify cells on Dirichlet bnd 
 Param.dof_f_dir =Grid.dof_f_xmax; % identify faces on Dirichlet bnd 
 Param.g = (Ta-th0)*(phi(Nx)*Cpi*rho_i+(1-phi(Nx))*Cpr*rho_r); 
 [B,N,fn] = build_bnd(Param,Grid,I,Nx); 
 
  H0=525*10^-9*0.357;%Al
  H1=0.063*150*10^-9;% Fe
  H2=29.17e-6*943e-9;  %40k
  H3=26.38e-6*44.3e-9;  %Th
  H4=568.7e-6* 6.27e-9;  %235U
  H5=94.65e-6*20.2e-9;  %238U
  lambda=log(2)/(0.717e6) ;%Al
  lambda1=log(2)/(1.5e6) ;
  lambda2=5.54*10^-10 ;  %40k
  lambda3=4.95*10^-11; 
  lambda4=9.85*10^-10; 
  lambda5=1.551*10^-10; 

  idw=zeros(Nx,1);
%phi: weight fraction
wa(1,:)=w;
phic(1,:)=phi;
Rj=0;
h10=(Tas-th0)*Cpi*rho_i; %the entholpy at melting temperature
%x=rho_r*phi./(rho_i*(1-phi)+rho_r*phi);
%Opropeties.x=x;
T=zeros(Nx,Lt);

d1 = @(w) (w==0);
dmix = @(w) 1-(w==0)-(w==1);
d2 = @(w) (w==1);

chi_w1=d1(w);
chi_m=dmix(w);
chi_w2=d2(w);
un=(1-phi)*Cpr*rho_r.*(To-th0)+phi*Cpi*rho_i.*(To-th0).*chi_w1...
    +phi.*(Cpi*rho_i*(Tas-th0)+Lh*w*rho_i).*chi_m...
    +phi.*(Cpi*(Tas-th0)*rho_i+Cpw.*(To-Tas)*rho_i+Lh*rho_w).*chi_w2;

uc=zeros(Nx,Lt);
uc(:,1)=un;


ic=1;
Tk=To;

%phi_x2=phi(1);
% phi_trgt=0.134156246097915;
phi_trgt=phis;

%for i=1:Lt0-1
i=0;
 while  phi(1)> phi_trgt  
    i=i+1
    
h1=h10*phi+(1-phi)*(Tas-th0)*Cpr*rho_r;
h2=h1+Lh*phi*rho_i; 
Cp=phi.*Cpi.*(1-w)*rho_i+phi.*Cpw.*w*rho_w+(1-phi)*Cpr*rho_r;


     Kappa=phi.*k_i.*(1-w)+phi.*k_w.*w+(1-phi).*k_r; %this Kappa=Conductivity
     pp=-1*ones(Nx+1,1);
      Kappam=comp_mean(Kappa,pp,Nx);
     %rho=phi*rho_i.*(1-w)+phi*rho_w.*w+(1-phi)*rho_r;
     %M=spdiags(I,0,Nx,Nx);
     dt=t1(i+1)-t1(i);
     difuT=(D*Kappam*G)*Tk;
    Fs0=un; 
    fs=rho_r.*(1-phi)*(H0/(lambda)*(exp(-lambda*t1(i))-exp(-lambda*t1(i+1)))...
      +H1/(lambda1)*(exp(-lambda1*t1(i))-exp(-lambda1*t1(i+1)))...
      +H2/(lambda2)*(exp(-lambda2*t1(i))-exp(-lambda2*t1(i+1)))...
      +H3/(lambda3)*(exp(-lambda3*t1(i))-exp(-lambda3*t1(i+1)))...
      +H4/(lambda4)*(exp(-lambda4*t1(i))-exp(-lambda4*t1(i+1)))...
      +H5/(lambda5)*(exp(-lambda5*t1(i))-exp(-lambda5*t1(i+1))));
     u = solve_lbvp(I,Fs0+fn+yr2s*(fs+dt*difuT),B,Param,N); 
  
  un=u;

  Rj=0;
      for j=1:Nx
        if un(j)<h2(j) && un(j)>h1(j)
         Tk(j)=Tas;
         w(j)=(un(j)-h1(j))/(h2(j)-h1(j));
         Rj=j;
         idw(j)=1;
        elseif un(j)<=h1(j)
        Tk(j)=th0+un(j)/Cp(j);
        w(j)=0;
         idw(j)=0;
        else
         Tk(j)=Tas+(un(j)-h2(j))/Cp(j);
         w(j)=1;
          %Rj=j;
           idw(j)=1;
        end
      end
        wold=w;
   %Rj;
   
%rearrange water
if  sum(idw)>=2 && w(1)<1
    nBotm= find(idw,1,'first');
    nTop= find(idw,1,'last');
     %0.134156246097915

    Vw=Grid.V'*(w.*phi);%volume of water we are transporting

    phiin=Grid.V(nBotm:nTop)'*phi(nBotm:nTop);  %volume of water and ice in target area

    V0=sum(Grid.V(nBotm:nTop));
    V_rock=V0-phiin;  %volume of rock phase
    phi_x2=(phiin-Vw)/(V0-Vw);  %volume fraction

if  phi_x2>=phi_trgt

    phi_x3=phi_x2;

    V_m=V_rock/(1-phi_x3);
  
 
         Vw_rs=0;

  %Vs=V0-Vw; %volume of solid in target
  w_rs=Vw_rs/(V_m-V_rock);
  Vp=Grid.V(nBotm:nTop);
  csVp=cumsum(Vp);
 %Rk=(Vs/(4/3*pi))^(1/3);
% Nj=floor(Rk/dx); % number of ice layers
  nVs=sum(V_m>csVp);
w(nBotm:nBotm+nVs-1)=w_rs;
phi(nBotm:nBotm+nVs-1)=phi_x3;
Vsrs=V_m-csVp(nVs); %residul solid volume
Vwk=Grid.V(nBotm+nVs)-Vsrs;
Vphik=Vwk+Vsrs*phi_x3;
wk=(Vwk+w_rs*phi_x3*Vsrs)/Vphik;
phik=Vphik/Grid.V(nBotm+nVs);
w(nBotm+nVs)=wk;
phi(nBotm+nVs)=phik;

if nTop-nBotm>=2
    w(nBotm+nVs+1:nTop)=1;
    phi(nBotm+nVs+1:nTop)=1;
end

chi_w1=d1(w);
chi_m=dmix(w);
chi_w2=d2(w);
%unold=un;
un=(1-phi)*Cpr*rho_r.*(Tk-th0)+phi*Cpi*rho_i.*(Tk-th0).*chi_w1...
    +phi.*(Cpi*rho_i*(Tas-th0)+Lh*w*rho_i).*chi_m...
    +phi.*(Cpi*(Tas-th0)*rho_i+Cpw.*(Tk-Tas)*rho_w+Lh*rho_w).*chi_w2;

else
    phi_x3=phi_trgt;
    V_m=V_rock/(1-phi_x3);
    Vw2=V0-V_m; %volume of water we actually transported

         Vw_rs=Vw-Vw2;   %liquid water volume left in the matrix

  %Vs=V0-Vw; %volume of solid in target
  w_rs=Vw_rs/(V_m-V_rock);
Vp=Grid.V(nBotm:nTop);
csVp=cumsum(Vp);
 %Rk=(Vs/(4/3*pi))^(1/3);
% Nj=floor(Rk/dx); % number of ice layers
  nVs=sum(V_m>csVp);
w(nBotm:nBotm+nVs-1)=w_rs;
phi(nBotm:nBotm+nVs-1)=phi_x3;
Vsrs=V_m-csVp(nVs); %residul solid volume
Vwk=Grid.V(nBotm+nVs)-Vsrs;
Vphik=Vwk+Vsrs*phi_x3;
wk=(Vwk+w_rs*phi_x3*Vsrs)/Vphik;
phik=Vphik/Grid.V(nBotm+nVs);
w(nBotm+nVs)=wk;
phi(nBotm+nVs)=phik;

if nTop-nBotm>=2
    w(nBotm+nVs+1:nTop)=1;
    phi(nBotm+nVs+1:nTop)=1;
end
 
chi_w1=d1(w);
chi_m=dmix(w);
chi_w2=d2(w);
unold1=un;
un=(1-phi)*Cpr*rho_r.*(Tk-th0)+phi*Cpi*rho_i.*(Tk-th0).*chi_w1...
    +phi.*(Cpi*rho_i*(Tas-th0)+Lh*w*rho_i).*chi_m...
    +phi.*(Cpi*(Tas-th0)*rho_i+Cpw.*(Tk-Tas)*rho_w+Lh*rho_w).*chi_w2;
unold2=un;



%%

end



end

 Opropeties.wa1=w;
 Opropeties.phi1=phi;
 Opropeties.T1=Tk;

      % write data
    if mod(i,Dt) == 0
        ic=ic+1;
        Flux.source(ic,:)=yr2s*fs; 
        uc(:,ic)=u; 
        T(:,ic)= Tk;
        wa(ic,:)=w;
        phic(ic,:)=phi;
        Kappac(ic,:)=Kappa;
        Cpc(ic,:)=Cp;
        ppc(ic,:)=pp';
        Flux.wf(ic,:)=Kappam*G*Tk;
        Flux.wF(ic,:)=yr2s*dt*(Kappam*G)*Tk;
    end
 end
 ti=i;
 T(:,1)=To;
 
 Opropeties.wa=wa;
 Opropeties.phi=phic;
 Opropeties.Kappac=Kappac;
 Opropeties.Cp=Cpc;   %Capacity
 Opropeties.pp=ppc; 
