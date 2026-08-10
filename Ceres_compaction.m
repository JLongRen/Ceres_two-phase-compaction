clear;
clc;
close all; 


ZKZ=1;  %to scale up or down solid viscosities
J=-0.5;   % solid viscosity mixture coeffcient
theta=0;  % theta for the numerical scheme 0:Forward Euler 1:Backward Euler
R=470*10^3;

Vo=4/3*pi*R^3;

T_sp=780.46; %K  pressure depended  P=100 MPa 

rho_i=1000;
rho_w=1000;
rho_sp=2545.87;  % density of serpentine from our thermodynamic method at 100 PA and 780K

Cpi=2100; %specific heat capacity for ice  J/(kg K)
Cpw=4200;  %specific heat capacity for water

Lh=333.55*1000;  %ice-water latent heat J/kg
Lhsp=3.15500E+05;  %latent heat of serpentinization under 100 MPa；  

CpSp=1201.65;

G_g=6.67*10^-11;
yr2s=365.25*24*3600;
ck=1e-15; 
Nx1=200;
Nt=180000;      %120000 for ck=1e-16;  180000 for ck=1e-15;  180000/3e5(mu-1) for ck=1e-14;  
% 300000 theta=1 for 1e-13; 9e5 mu1 15e5 mu-1 theta=1   for k=1e-12

detN=300; % detN about equal to Nt/600
tf=3e6;   % accretion time (yr)
te=5000e6;  % simulation end time (yr)
YrperN=te/Nt;
rho_m0=2162;
phiTw=0.25*rho_m0/rho_w; %volume fraction of total water assume no volume change at reactions
rho_olx=(rho_m0-rho_w*phiTw)/(1-phiTw);
rho_ol=rho_olx;

wsp_const=(rho_ol-rho_sp)/(rho_ol-rho_w);  %water volume fraction in serpenite assuming DV=0
phiv=(phiTw-wsp_const)/(1-wsp_const);

wsp_C=wsp_const*rho_w/rho_sp;
wsp_isp=rho_sp/(1-wsp_const); 
nx=3;

muf=8e-4;  % fluid viscosity for the characteristic values calculation
Zeta0=1e20;
m=-1;
phi(1:Nx1)=phiv;
phi=phi';
ZetaC=Zeta0*phiv^-m; %characteristic solid viscosity
delta=sqrt(ck*ZetaC/muf);
drho=rho_w-rho_sp;
pc=-drho*delta;

 tc=ZetaC/pc;   %compaction time scale  s
 tcy=tc/yr2s;   %compaction time scale  yrs
 vcr=delta/tc;
 qrc=ck*pc/(muf*delta);
 %%  thermal properties
paraphysic.k_i=2.3;
paraphysic.k_sp=2.11; %1
paraphysic.k_ol=2.56;
paraphysic.k_w=0.6;
paraphysic.Cpr=CpSp;
paraphysic.Cpi=2100;
paraphysic.Cpw=4200;
paraphysic.Lh=Lh;

paraphysic.rho_i=1000;
paraphysic.rho_w=1000;
paraphysic.rho_r=rho_sp;

%%  Grid (dimensionless) defination
Grid.xmin = 0; Grid.xmax =R/delta; Grid.Nx =Nx1; %
Grid.geo=2; 
Grid = build_grid(Grid); 
dx0=Grid.dx;
grx0=Grid.xc;
grxf=Grid.xf;
V_ra=4/3*pi*((grx0+0.5*dx0).^3-(grx0-0.5*dx0).^3)';

[Dr,Gr,Ir]=build_ops(Grid);
%% accretional temperaute profile
G=6.67*10^-11;
 To=155; %K
ha=0.2;
theta2=4;
Cp_acc=1178;  %J/(kg*K)
T_acc=@(r)  ha/Cp_acc*(4/3*pi*rho_m0*G*r.^2)*(1+0.5*theta2);
Tacc1=T_acc(grx0*delta)'+To;
Tacc1(Nx1)=To;
%%  initial accumulation
tw=4e6;
Nt0=12000;
tor0=linspace(tf,tw,Nt0);
phis=0.15;
[ur0,T0,proper0,Flux0,ti]=wout5(R,Tacc1,tor0,Nx1,phi,To,0*ones(Nx1,1),10,paraphysic,phis);
phiIX=proper0.phi1;
tf2=tor0(ti+1);

t1=linspace(tf2,te,Nt+1);

%t2=union(t1(1:detN:Nt+1),t1(1:2:2*detN)); %to restore the early stage in detail
t2=t1(1:detN:Nt+1);
[Lia,Locb] = ismember(t2,t1);
Ntx2=size(t2,2);

dt=yr2s*(te-tf2)/Nt/tc; 

phi=phiIX;
X_phi=(phi==1);
nb_ocean= find(X_phi,1,'first'); 
topo=find(X_phi,1,'last'); 

rho_ra0=phi(1)*rho_w+(1-phi(1))*rho_sp;
m0=4/3*pi*rho_ra0*Grid.xmin^3;
V_raC=4/3*pi*((grx0).^3-([Grid.xmin, grx0(1:Nx1-1)].^3))';

rho_mean=(phi.*rho_w+(1-phi).*rho_sp)'*V_ra/sum(V_ra); %check the mean density
drho=rho_w-rho_sp;

vf_dgh=ones(Nx1,1);   %serp(hydrated) degree from 0 to 1
vf_aw=vf_dgh*wsp_const.*(1-phi)+phi;  %initial free+combined water volume fraction


%%  w: water mass fraction /(ice+water)
th0=155; %K surface T
Tas=273.15;

d1 = @(w) (w==0);
dmix = @(w) 1-(w==0)-(w==1);
d2 = @(w) (w==1);

w1=zeros(Nx1,1);    %initially all ice
w1(1:topo)=1;
un=@(T,phi,w) ((1-phi)*rho_sp*CpSp).*(T-th0)+phi*rho_i*Cpi.*(T-th0).*(1-(w==1))...
             +(Lh*rho_i*w.*phi).*(1-(w==0))...
             +(phi.*rho_i*Cpi*(Tas-th0)+phi.*rho_w*Cpw.*(T-Tas)).*(w==1);
un1=un(300,phi,w1);

h10=(Tas-th0)*Cpi*rho_i; %the enthalpy at melting temperature
hw3=h10+Lh*rho_i+(T_sp-Tas)*Cpw*rho_w;  %entholpy of unit voluem water at T_sp
hsp3=(T_sp-th0)*CpSp*rho_sp;
hol3=(hsp3+Lhsp*rho_sp-wsp_const*hw3)/(1-wsp_const);  %enthalpy of olivine per unit volume*K



CpOl=hol3/(rho_ol*(T_sp-th0));
%%  fluid viscosity
mufVecot = readmatrix('D:\Projects\Ceres compaction\compaction\mufdata'...
    ,'Sheet','100MPa');  %change it to the direction of the fluid viscosity data 
ln_muf=log(mufVecot);
Ti=linspace(0.1,700,100)+Tas;  % in K

pply = polyfit(Ti,ln_muf,4);
%% solid viscosity
Tus=th0; %155  %K
dgrain=20; % 1e-6  m  
Zeta_sp=1e19;  %Pa*s
Ex=1e-15;
ZetadT=@(T_Du) (10^5.62*10^-6*dgrain^(-3).*exp(-240e3./(8.31*T_Du))).^(-1);

Zeta_Q=@(T,n,E) 0.5*E^((1-n)/n)*3^(-(n+1)/(2*n))*(exp(13.6)/2)^(-1/n)*1e6...
        .*exp(480e3./(n*8.31*T));                   % Miller 2012

Zeta_spT=@(T,n,E) E^((1-n)/n)*(exp(86/n)).*exp((8.9e3+3.2*100)./(n*8.31*T));  %n=3.8  % Hilairet 2007

tcpdt=@(TT) -1/drho*sqrt(muf*ZetadT(TT)./ck)/(tc*dt); %haven't consider phi dependence 

deT=500+Tas-Tus;
nTX=100;

T_SI(1:Nx1)=Tas;
T_SI(Nx1)=Tus;

ZeK=ZetadT(T_SI);
ZeKo=phi.^m.*ZeK';
tcR=tcpdt(T_SI);
%% define the recording matrices 
phi_c=zeros(Nx1,Ntx2);
v_sp_c=zeros(Nx1,Ntx2);  %volume fraction of srp
v_ol_c=zeros(Nx1,Ntx2);
p_c=zeros(Nx1,Ntx2);
pf_c=zeros(Nx1,Ntx2);
v_c=zeros(Nx1+1,Ntx2);
v_ec=zeros(Nx1,Ntx2);
qr_c=zeros(Nx1+1,Ntx2);
H_c=zeros(Nx1,Ntx2);
C_c=zeros(Nx1,Ntx2);
T_c=zeros(Nx1,Ntx2);
w_c=zeros(Nx1,Ntx2);
ps_c=zeros(Nx1,Ntx2);
vf_sp_c=zeros(Nx1,Ntx2);

Zeta_c=zeros(Nx1,Ntx2);
delta_c=zeros(Nx1,Ntx2);
Zeta_c(:,1)=ZeKo;
w_c(:,1)=w1;
vf_sp_c(:,1)=1;
v_ol_c(:,1)=0;
deVN=zeros(Ntx2,1);
deWN=zeros(Ntx2,1);
deW1N=zeros(Ntx2,1);
DWN=zeros(Nx1,Ntx2);
NNzone=zeros(Ntx2,1);
NNzone(1)=1;

p_s=zeros(Nx1,Ntx2);
p_f=zeros(Nx1,Ntx2);

V_sys_c=zeros(Nx1,Ntx2);
LL_c=zeros(Nx1,Ntx2);
vf_aw_c=zeros(Nx1,Ntx2);
Ol_ex_c=zeros(Nx1,Ntx2);
w_ex_c=zeros(Nx1,Ntx2);

stainR_c=zeros(Nx1,Ntx2);

Tc=deT; 
nzN=zeros(Ntx2,1);
nzN2=zeros(Ntx2,1);
T=(T_SI-Tus)'/Tc;  %dimensionless T
Ho=un(Tc,phiv,1);
H=un(T_SI',phi,w1)/Ho;
Hm=H;

h1n=h10/Ho;
h2n=h1n+(Lh*rho_i)/Ho; 

H_c(:,1)=H;
T_c(:,1)=(T_SI-Tus)/Tc;
phi_c(:,1)=phi;

v_sfaceh=zeros(Nx1+1,1);
grav=zeros(Nx1+1,1);
%% moment of inertia
M=4/3*pi*rho_m0*R^3; %kg
grxDm=grx0*delta;
rho_x=phi.*rho_w+(1-phi).*rho_sp;
Isph=8/3*pi*(grxDm.^4*rho_x*dx0)*delta;
Irf=2/5*M*R^2;
IoMR=Isph/(M*R^2);
Iorf=Isph/Irf;

  %%  Radioactive heating constants
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
  %%
nzN(1)=Nx1;
nzN2(1)=nb_ocean;
deV=0;
p=zeros(Nx1,1);

S_Phi0=sum(delta^3*V_ra'*phi/Vo,1);
vf_olnw=zeros(Nx1,1);
vf_ol=zeros(Nx1,1);  % volume fraction of ol 
V_Sp=(1-phi).*ones(Nx1,1);  % volume fraction of sp 

C=vf_aw;   %dimensionless total water composition 
C_c(:,1)=C;

spR_dgr=zeros(Nx1,1); %deserpen degree
Ol_ex=zeros(Nx1,1);

v_sp1=vf_aw/wsp_const;   %volume fraction of sp when all free water is consumed
v_sp2=(1-vf_aw)/(1-wsp_const);   %volume fraction of sp when all olivine is consumed

v_spTh=min(v_sp1,v_sp2);   % potential max volume fraction of sp
v_sp=v_spTh;
v_ol=1-v_sp-phi;
v_sp_c(:,1)=v_sp;

h1=(h10*(1-v_spTh)+v_spTh*(Tas-th0)*CpSp*rho_sp)/Ho;   %left end at Tmelt
h2=h1+(Lh*(1-v_spTh)*rho_i)/Ho; %Right end at Tmelt
h3=1/Ho*(v_spTh*CpSp*rho_sp*(T_sp-th0)+(1-v_spTh).*hw3.*(vf_aw>=wsp_const)+(1-v_spTh).*hol3.*(vf_aw<wsp_const));  %left end at T_sp
h4=1/Ho*(hw3*vf_aw+hol3*(1-vf_aw));  %Right end at T_sp

Cp=(phi.*Cpi.*(1-w1)*rho_i+phi.*Cpw.*w1*rho_w+V_Sp*CpSp*rho_sp+vf_ol*CpOl*rho_ol)/Ho;
phi2=phi;
vf_aw_c(:,1)=vf_aw;
Vc=sum(Grid.V);
totW0=Grid.V'*vf_aw/Vc;

V_uphalf=4/3*pi*(grxf(2:end).^3-grx0.^3)';
V_lowhalf=4/3*pi*(grx0.^3-grxf(1:end-1).^3)';
cnwt=zeros(Nx1,1);
psD=zeros(Nx1,1);
ps=zeros(Nx1,1);
D_C=zeros(Nx1,1);

%%  main calculation
j=1;

for i=1:Nt

LocbX=Locb(j+1);
v_sfaceh=zeros(Nx1+1,1);  %update here to avoid use faulse vs & qr in thermal calculation
qrh=zeros(Nx1+1,1);
pfull=zeros(Nx1+1,1);

D_C=wsp_const*(T<=T_sp); 
nndc=zeros(Nx1,1);
   nndc(phi<1)=1;
   nndc=logical(nndc);

D_C(nndc)=wsp_const*rho_w*v_sp(nndc)./(1-phi(nndc))/rho_w ;

D_H=CpSp*rho_sp/(Cpw*rho_w)*ones(Nx1,1);  
nndh=zeros(Nx1,1);
   nndh(phi<1)=1;
   nndh=logical(nndh);
D_H(nndh)=(CpSp*rho_sp*v_sp(nndh)+CpOl*rho_ol*v_ol(nndh))./(1-phi(nndh))./(Cpw*rho_w) ;

    Hold=H;
    Told=T;
    phiold=phi;
    Cold=C;
    phi2old=phi2;
    v_spold=v_sp;
    v_olold=v_ol;
    totWold=totW0;
    vf_awold=vf_aw;

phiIR=phi.*w1;
Tdm=T'*Tc+Tus;
tcR=phi'.^(0.5*(m-nx)).*tcpdt(Tdm);

S_phi=(phiIR>=1-1e-6);
X_phi=(phi==1);
S_W=(w1<=1e-6);
nb_ocean= find(X_phi,1,'first'); 
n1phi= find(S_phi,1,'first');  % pure liquid water ocean
n2phi= find(S_W,1,'first');
zN2= find(S_phi,1,'last');

zN=min(n2phi,nb_ocean)-1; %the upper bound of compaction zone(s)

% count how many disconnected compaction zones
phi_ny=phi(1:zN);
x1=zeros(zN,1);
x1(phi_ny>0)=1;
xm=[0; x1];
sqrX=(xm(2:end)-xm(1:end-1)).^2;

njump=sum(sqrX);
xx2=linspace(1,zN,zN)';
yy2=xx2.*sqrX;
yy2(yy2==0) = [];

N_zone=ceil(njump/2); % the number of compaction zones
i_bc=1:njump;
i_left=yy2(mod(i_bc,2)==1);

if mod(njump,2)==1
        i_right=[yy2(mod(i_bc,2)==0)-1; zN];
else
        i_right=yy2(mod(i_bc,2)==0)-1;
end

nxx=zeros(Nx1,1);

% update the gravity field
rho_ra=phi*rho_w+v_sp*rho_sp+v_ol*rho_ol;
m_ra=V_ra.*rho_ra;
acm_ra=m0+cumsum(m_ra);
acm_raC=acm_ra-rho_ra.*4/3*pi.*(grxf(2:end)'.^3-grx0'.^3);
grav(2:Nx1+1)=delta*G_g*acm_ra./(grxf(2:Nx1+1)'.^2);  %re-dimensional for g
grav(1)=0;
gravC=delta*G_g*acm_raC./(grx0'.^2);

v_sp1=[v_sp(2:end);0];  %don't need the center
v_sp2=(v_sp1.*V_lowhalf+v_sp.*V_uphalf)./V_ra;

v_ol1=[v_ol(2:end);0];  %don't need the center
v_ol2=(v_ol1.*V_lowhalf+v_ol.*V_uphalf)./V_ra;

cnwt0=(v_sp2*rho_sp+v_ol2*rho_ol).*gravC*delta;
cnwt(1:nb_ocean-1)=cnwt0(1:nb_ocean-1);
cnwt(nb_ocean:Nx1)=cnwt0(nb_ocean:Nx1)+(phi(nb_ocean:Nx1)*rho_w).*gravC(nb_ocean:Nx1)*delta;
psf=cumsum(flip(cnwt));    %mass of solid phase per area, dimensional
psD0=flip(psf);  

psD=psD0;
ps=psD/pc;

% if i==1
%     ps_c(:,1)=psD;
% end

% if zN>0
% phi(zN)=min(phi(zN),0.6);
% end

for  i_1=1:N_zone
    
    zN_L=i_left(i_1);
    zN_R=i_right(i_1);
   nxx(zN_L:zN_R)=1;
    zN_w=[zN_L zN_R];
phi_nz=phi(i_left(i_1):i_right(i_1));
nLL=zN_R+1-zN_L;
 
xc=grx0(zN_L:zN_R);
xf=grxf(zN_L:zN_R+1);

Da=Dr(zN_L:zN_R,zN_L:zN_R+1);
Ga=Gr(zN_L:zN_R+1,zN_L:zN_R);
Ga(1,:)=0;
Ga(nLL+1,:)=0;
Ia=Ir(zN_L:zN_R,zN_L:zN_R);


mu_T= exp(polyval(pply,Tdm)); %for fluid viscosity

Zeta_ol=Zeta_Q(Tdm,3.5,Ex);

Zeta_spX=Zeta_spT(Tdm,3.8,Ex);
vv_ol=v_ol'./(1-phi');
vv_sp=v_sp'./(1-phi');
Zeta_mix=(vv_ol.*Zeta_ol.^J+vv_sp.*Zeta_spX.^J).^(1/J); %except J=0

ZeK=Zeta_mix./mu_T;
AA=ZeK(zN_L:zN_R)'./(ZetaC/muf);
Zeta=ZKZ*AA.*phi_nz.^m;

k=phi_nz.^nx;
K_0=k;
K=comp_mean(K_0,-1*ones(nLL+1,1),nLL);
Isan=spdiags(1./Zeta,0,nLL,nLL);
     L = -Da*K*Ga+Isan;
     RP=-Da*(K*grav(zN_L:zN_R+1));

        Param.dof_neu =zN_L ;
        Param.dof_f_neu =zN_L; 
        Param.qb=-grav(Param.dof_f_neu).*k(1)/(mu_T(zN_L)/muf);

        Param.dof_dir=zN_R ;
        Param.dof_f_dir=zN_R+1; 
        Param.g=0;

 [Bp,Np,fn_p] = build_bndZ(Param,Grid,Ia,zN_L,zN_R);
 p= solve_lbvp(L,fn_p+RP,Bp,Param,Np); 
  
 pfull(zN_L:zN_R)=p;
flux = @(u) -K*(Ga*u-grav(zN_L:zN_R+1));
res = @(u,cell) L(cell,:)*u - RP(cell);  

 qr = comp_flux_genY(flux,res,p,Grid,Param,zN_L,zN_R);

pR=p*pc;
qr(1)=0;
qr(nLL+1)=0;
%--------------  
 Paramu.dof_neu =[] ;     %compaction length delta is inside of uc us=-drho*delta*k0/mi\uf
 Paramu.dof_f_neu =[]; 
 Paramu.qb=[];

      Lu = -Da*Ga;
      Ru= p./Zeta;

 Paramu.dof_dir =zN_R ; % identify cells on Dirichlet bnd 
 Paramu.dof_f_dir =zN_R+1; % identify faces on Dirichlet bnd
 Paramu.g = 30; 
   [Bu,Nu,fn] = build_bndZ(Paramu,Grid,Ia,zN_L,zN_R);
  u = solve_lbvp(Lu,Ru+fn,Bu,Paramu,Nu); %velocity potential

flux =  @(u) -Ga*u;
res = @(u,cell) Lu(cell,:)*u - Ru(cell);
vD = comp_flux_genY(flux,res,u,Grid,Paramu,zN_L,zN_R);
v_sface=-qr;
v_sfaceh(zN_L:zN_R+1)=v_sface;
stainR=-qrc/delta*Dr*v_sfaceh;
qrh(zN_L:zN_R+1)=qr;

%--------------

      ParamC.dof_dir=[];% identify cells on Dirichlet bnd 
      ParamC.dof_f_dir =[]; % identify faces on Dirichlet bnd 
      ParamC.g=[];

      ParamC.dof_neu =[] ;
      ParamC.dof_f_neu =[]; 
      ParamC.qb=[];

  A1 = flux_upwind(qrh(zN_L:zN_R+1),nLL,nLL+1)*spdiags(1./(phi(zN_L:zN_R)+(1-phi(zN_L:zN_R)).*D_C(zN_L:zN_R)),0,nLL,nLL);
  A2 = flux_upwind(v_sfaceh(zN_L:zN_R+1),nLL,nLL+1);

  ve_mC=A1+A2;
  [Lim,Lex] = build_ade_opsA(Da,ve_mC,Ia,dt) ;
   
  [BC,NC,fnC] = build_bndZ(ParamC,Grid,Ia,zN_L,zN_R);
  IM=Lim(theta); 
  EX=Lex(theta); 

  C(zN_L:zN_R) = solve_lbvp(IM,EX*C(zN_L:zN_R),BC,ParamC,NC);
%---------------------------------------------------

    Paramphi.dof_dir=[] ;
    Paramphi.dof_f_dir=[]; 
    Paramphi.g=[];

    Paramphi.dof_neu =[] ;
    Paramphi.dof_f_neu =[]; 
    Paramphi.qb=[];

    [Bphi,Nphi,fnphi] = build_bndZ(Paramphi,Grid,Ia,zN_L,zN_R);
    [Lim,Lex] = build_ade_ops(nLL,nLL+1,Da,v_sfaceh,Ia,dt) ;

    Dvs=Da*v_sface;
    IM=Lim(theta);  
    EX=Lex(theta);
    Fbnd=zeros(nLL,1);
    Fs0=EX*phi(zN_L:zN_R)+dt*Dvs+Fbnd; 
    
    phi2(zN_L:zN_R) = solve_lbvp(IM,Fs0+fnphi,Bphi,Paramphi,Nphi);
end

if N_zone>0
%both phi and C at top need to be evaluated seperately
nxx=logical(nxx);

phi(nxx)=phi2(nxx);
C_tmp=C;

Dphi=phi-phiold;
deV=-sum(Dphi(1:zN).*Grid.V(1:zN));  %liquid volume that moves upward from compaction
vf_aw=C;
vf_awcmp=vf_aw;
DW0=vf_aw-vf_awold;
deW=-sum(DW0(1:zN-1).*Grid.V(1:zN-1));

deW_N=zeros(N_zone,1);
deV_N=zeros(N_zone,1);
for i_1=1:N_zone
    zL=i_left(i_1);
    zR=i_right(i_1);
    deW_N(i_1)=-sum(DW0(zL:zR-1).*Grid.V(zL:zR-1));
    deV_N(i_1)=-sum(Dphi(zL:zR-1).*Grid.V(zL:zR-1));
    
    ndeW_sum=[flip(cumsum(DW0(zL:zR-1).*Grid.V(zL:zR-1))); 0];

    % fill the compaction top with water:
    ndeV_sum=[flip(cumsum(-Dphi(zL:zR-1).*Grid.V(zL:zR-1))); 0]; 
    
    vf1_zN=vf_awold(zR)*Grid.V(zR); % how much more free water could zN cell give
    vf2_zN=(1-vf_awold(zR))*Grid.V(zR); % how much more free water could zN cell take, dimensionless
    phi1_zN=phiold(zR)*Grid.V(zR);  
    phi2_zN=(1-phiold(zR))*Grid.V(zR);

    phi_sum=cumsum(flip(phiold(zL:zR)).*flip(Grid.V(zL:zR)));
    vwbacksp=phiold(zL:zR).*wsp_const.*v_spold(zL:zR)./(v_spold(zL:zR)+v_olold(zL:zR));
    phi_mod=phiold(zL:zR)-vwbacksp;
    phi_msum=cumsum(flip(phi_mod).*flip(Grid.V(zL:zR)));
    phiREAL_sum=cumsum(flip(phiold(zL:zR)).*flip(Grid.V(zL:zR)));
    newzN=find((phi_msum>=ndeW_sum),1,'first');
    
    if isnan(newzN)
            msg = 'wrong deW';
            error(msg)
    end
    
    if deW_N(i_1)>vf2_zN
        deWnext=deW_N(i_1)-vf2_zN;
        C(zR)=1;
        vfx=deWnext/Grid.V(zR-1)+vf_aw(zR-1);
        C(zR-1)=vfx;
    else 
        vfx=deW_N(i_1)/Grid.V(zR)+vf_awold(zR);
        C(zR)=vfx;
    end

    if deV_N(i_1)>phi2_zN
        deVnext=deV_N(i_1)-phi2_zN;
        phi(zR)=1;
        phi(zR-1)=phi(zR-1)+deVnext/Grid.V(zR-1);
    else
        phi(zR)=phiold(zR)+deV_N(i_1)/Grid.V(zR);
    end

    if  newzN>1
        drVnext=ndeW_sum(newzN)-phi_msum(newzN-1); 
        phi(zR-newzN+2:zR)=0;
        vfx0=wsp_const*v_spold(zR-newzN+2:zR)./(v_spold(zR-newzN+2:zR)+v_olold(zR-newzN+2:zR));
        C(zR-newzN+2:zR)=vfx0;
        deVX=-sum(Dphi(zL:zR-newzN).*Grid.V(zL:zR-newzN));% a potential bug
        phi(zR-newzN+1)=phiold(zR-newzN+1)+(phiREAL_sum(newzN-1)+deVX)./Grid.V(zR-newzN+1);
        vfx1=-drVnext/Grid.V(zR-newzN+1)+vf_awold(zR-newzN+1);
        C(zR-newzN+1)=vfx1;
    end
end


vf_aw=C;  %
v_sp(nxx)=(vf_aw(nxx)-phi(nxx))/wsp_const;  %volume fraction of all water component
v_ol(nxx)=1-v_sp(nxx)-phi(nxx);

DW1=vf_aw-vf_awold;
deW1=-sum(DW1(1:zN).*Grid.V(1:zN));

end

Cp=(phi.*Cpi.*(1-w1)*rho_i+phi.*Cpw.*w1*rho_w+v_sp*CpSp*rho_sp+v_ol*CpOl*rho_ol)/Ho;   %!!
%---the compaction ends here----------------------------------
if   max(vf_aw)>1 
    msg = 'too much water';
    error(msg)
end

if   min(vf_aw)<0 
    msg = 'too less water';
    error(msg)
end

if   max(phi)>1 || min(phi)<0 
    msg = 'wrong phi';
    error(msg)
end

nn=vf_aw<1;
nW=vf_aw==1;

phi0=(vf_aw-wsp_const*(1-vf_aw)./(1-wsp_const)).*(vf_aw>wsp_const); %porosity below h3

v_sp1=vf_aw/wsp_const;   %volume fraction of sp when all free water is consumed
v_sp2=(1-vf_aw)/(1-wsp_const);   %volume fraction of sp when all olivine is consumed
v_spTh=min(v_sp1,v_sp2);   % potential max volume fraction of sp

h1=(h10*(1-v_spTh)+v_spTh*(Tas-th0)*CpSp*rho_sp)/Ho;   %left end at Tmelt
h2=h1+(Lh*(1-v_spTh)*rho_i)/Ho; %Right end at Tmelt
h3=1/Ho*(v_spTh*CpSp*rho_sp*(T_sp-th0)+(1-v_spTh).*hw3.*(vf_aw>wsp_const)+(1-v_spTh).*hol3.*(vf_aw<=wsp_const));  %left end at T_sp
h4=1/Ho*(hw3*vf_aw+hol3*(1-vf_aw));  %Right end at T_sp
N_ocean=(phiold>=1-1e-6).*(w1==1);  %pure liquid water ocean
 %----------------------------------------------
 %ocean layer  assume ocean is well mixed (Tmelt for all)
sumO=sum(N_ocean);
T(N_ocean==1)=(Tas-Tus)/Tc;
H(N_ocean==1)=un(Tas,1,1)/Ho; 

 %----------------------------------------------------
if sum(N_ocean)>0
    yNa=n1phi;  Tub=Tas;   

    Param3.dof_dir=[];% identify cells on Dirichlet bnd 
    Param3.dof_f_dir =[]; % identify faces on Dirichlet bnd
    Param3.g=[];

    Param3.dof_neu =[] ;
    Param3.dof_f_neu =[]; 
    Param3.qb=[];

else
    yNa=Nx1;  Tub=Tus;

      Param3.dof_dir=yNa;% identify cells on Dirichlet bnd 
      Param3.dof_f_dir =yNa+1; % identify faces on Dirichlet bnd 
      Param3.g=un(Tub,phiold(yNa),w1(yNa))/Ho;

      Param3.dof_neu =[] ;
      Param3.dof_f_neu =[]; 
      Param3.qb=[];

end
Day=Dr(1:yNa,1:yNa+1);
Gay=Gr(1:yNa+1,1:yNa);
Gay(yNa+1,:)=0;
Iay=Ir(1:yNa,1:yNa);
 
  clear phimid
  phimid=0.5*(phi(1:yNa)+phiold(1:yNa));
  A1 = flux_upwind(qrh(1:yNa+1),yNa,yNa+1)*spdiags(1./(phimid+(1-phimid).*D_H(1:yNa)),0,yNa,yNa);
  A2 = flux_upwind(v_sfaceh(1:yNa+1),yNa,yNa+1);
  ve_m=A1+A2;
  ve_v=diag(ve_m,-1)+diag(ve_m);  
  [Lim,Lex] = build_ade_opsA(Day,ve_m,Iay,dt) ;
  pp=-1*ones(yNa+1,1);

  Kappawr=w1(1:yNa).*phimid(1:yNa).*paraphysic.k_w+(1-w1(1:yNa)).*paraphysic.k_i...
      +v_ol(1:yNa).*paraphysic.k_ol+v_sp(1:yNa).*paraphysic.k_sp;     
  Kappam=comp_mean(Kappawr,pp,yNa);
   
     [B3,N3,fn3] = build_bnd(Param3,Grid,Iay,yNa);
     difuT=(dt*tc/delta^2)*Tc*(Day*Kappam*Gay)*T(1:yNa)/Ho;
        
      fsz=(v_sp(1:yNa)*rho_sp+v_ol(1:yNa)*wsp_isp)*(H0/(lambda)*(exp(-lambda*t1(i))-exp(-lambda*t1(i+1)))...
            +H1/(lambda1)*(exp(-lambda1*t1(i))-exp(-lambda1*t1(i+1)))...
            +H2/(lambda2)*(exp(-lambda2*t1(i))-exp(-lambda2*t1(i+1)))...
            +H3/(lambda3)*(exp(-lambda3*t1(i))-exp(-lambda3*t1(i+1)))...
            +H4/(lambda4)*(exp(-lambda4*t1(i))-exp(-lambda4*t1(i+1)))...
            +H5/(lambda5)*(exp(-lambda5*t1(i))-exp(-lambda5*t1(i+1))));

      fs=tc/Ho*fsz/(t1(i+1)-t1(i));
      IM=Lim(theta); 
      EX=Lex(theta); 
      Fs=EX*H(1:yNa)+difuT+dt*fs; 

  Hm= solve_lbvp(IM,Fs,B3,Param3,N3); 
  H(1:yNa)=Hm; 

 
 T_dmL=(vf_aw(1:yNa)<=wsp_const).* ((Hm<=h3(1:yNa)).*(th0+Hm./Cp(1:yNa)) + T_sp.*(Hm>h3(1:yNa)).*(Hm<=h4(1:yNa))...
     +(T_sp+(Hm-h4(1:yNa))./Cp(1:yNa)).*(Hm>h4(1:yNa)));     

 T_dmR=(vf_aw(1:yNa)>wsp_const).* ((Hm<=h1(1:yNa)).*(th0+Hm./Cp(1:yNa)) + Tas.*(Hm>h1(1:yNa)).*(Hm<=h2(1:yNa))... 
     +(Tas+(Hm-h2(1:yNa))./Cp(1:yNa)).*(Hm>h2(1:yNa)).*(Hm<=h3(1:yNa)) + T_sp.*(Hm>h3(1:yNa)).*(Hm<=h4(1:yNa))...
     +(T_sp+(Hm-h4(1:yNa))./Cp(1:yNa)).*(Hm>h4(1:yNa))); 

 T_dm=T_dmL+T_dmR;
 Tm=(T_dm-Tus)/Tc;   % different reference T are used for Td and Hd

 flux =  @(u) -Tc/delta*Kappam*Gay*u;
 res = @(u,cell) Iay(cell,:)*u - Fs(cell);
 Qfluxm = comp_flux_genzq(flux,res,Tm,Hm,Grid,Param3,yNa,0);

 Qz0=Qfluxm(yNa+1);  %=0 in this case
 Qz1=@(n) Grid.A(yNa+1)*Qz0/Grid.A(n);

%------------------------------------------------------------
if sum(N_ocean)>0  %if there's no ocean layer, than ice and trans layers are included in the main part above

       N_ishell=(w1(zN2+1:end)==0);
       N_inter=(1-(w1==0)-(w1==1));

        szNic=sum(N_ishell);
        icN=Nx1-szNic+1;
        Daic=Dr(icN:Nx1,icN:Nx1+1);
        Gaic=Gr(icN:Nx1+1,icN:Nx1);
        Gaic(1,:)=0;
        Iaic=Ir(icN:Nx1,icN:Nx1);
        Hic=H(icN:Nx1);

        Paramic.dof_dir=[icN,Nx1]';% identify cells on Dirichlet bnd 
        Paramic.dof_f_dir =[icN,Nx1+1]'; % identify faces on Dirichlet bnd
        Paramic.g=[un(Tas,1,0)/Ho, un(Tus,1,0)/Ho]';
        Paramic.dof_neu =[] ;
        Paramic.dof_f_neu =[]; 
        Paramic.qb=[];
        [Bic,Nic,fnic] = build_bndq(Paramic,Grid,Iaic,szNic,icN-1);

        Kappa=phiold(icN:Nx1)*paraphysic.k_i+(1-phiold(icN:Nx1)).*paraphysic.k_sp; %this Kappa=Conductivity

        pp=-1*ones(szNic+1,1);
        Kappaic=comp_mean(Kappa,pp,szNic);
     
        difuT=(dt*tc/delta^2)*Tc*(Daic*Kappaic*Gaic)*T(icN:Nx1)/Ho;
     
        fscr=(v_sp(icN:Nx1)*rho_sp+v_ol(icN:Nx1)*wsp_isp)*(H0/(lambda)*(exp(-lambda*t1(i))-exp(-lambda*t1(i+1)))...
             +H1/(lambda1)*(exp(-lambda1*t1(i))-exp(-lambda1*t1(i+1)))...
             +H2/(lambda2)*(exp(-lambda2*t1(i))-exp(-lambda2*t1(i+1)))...
             +H3/(lambda3)*(exp(-lambda3*t1(i))-exp(-lambda3*t1(i+1)))...
             +H4/(lambda4)*(exp(-lambda4*t1(i))-exp(-lambda4*t1(i+1)))...
             +H5/(lambda5)*(exp(-lambda5*t1(i))-exp(-lambda5*t1(i+1))));

        fs=tc/Ho*fscr/(t1(i+1)-t1(i)); 
        Fs=Iaic*Hic+difuT+dt*fs; 
    
        Hic = solve_lbvp(Iaic,Fs,Bic,Paramic,Nic); 
        Tic=Ho*Hic./(rho_i*Tc*Cpi);
        flux =  @(u) -Tc/delta*Kappaic*Gaic*u;
        res = @(u,cell) Iaic(cell,:)*u - Fs(cell);
        Qflux = comp_flux_genzq(flux,res,Tic,Hic,Grid,Paramic,szNic,icN-1);
        H(icN:Nx1)=Hic; 

        szNin=sum(N_inter);
        inN=Nx1-szNic-szNin+1;  %lower cell
        dH=Grid.A(yNa+1)*Qz0-Qflux(1)*Grid.A(icN); 

if  szNin>0   

    nCH=icN-1;
    Holdx=  H(icN-1);
    H(nCH)=Holdx+dH/Grid.V(nCH);
    hx1=h1n-H(nCH); 
    hx2=H(nCH)-h2n;
    
    if  hx1>0   %freeze icN-2
      
        H(icN-1)=h1n;      
        H(icN-2)=H(icN-2)-hx1*Grid.V(icN-2+1)/Grid.V(icN-2);

    end

    if  hx2>0  %melt icN
        H(icN-1)=h2n;
        H(icN)=H(icN)+hx2*Grid.V(icN-1)/Grid.V(icN);
    end

else   % for ocean and iceshell directly contact  freeze icN-1 or melt icN
    if   dH>0  %melt icN
        nCH=icN;
    else  %freeze icN-1
        nCH=icN-1;
    end 
        H(nCH)=H(nCH)+dH/Grid.V(nCH);
   
end
end

vf_dghold=vf_dgh;


spR_dgr(nn)=((H(nn)-h3(nn))./(h4(nn)-h3(nn))).*(H(nn)<=h4(nn)).*(H(nn)>h3(nn))+(H(nn)>h4(nn));

vf_dgh(nn)=1-spR_dgr(nn);

phi(nn)=phi0(nn).*(vf_aw(nn)>wsp_const).*(H(nn)<=h3(nn))+(phi0(nn)+(vf_aw(nn)-phi0(nn)).*spR_dgr(nn)).*(H(nn)<=h4(nn)).*(H(nn)>h3(nn))...
 +vf_aw(nn).*(H(nn)>h4(nn));

phi(nW)=1;

Ol_ex(nn)=(vf_aw(nn)<=wsp_const).*(1-vf_aw(nn)/wsp_const);
v_ol(nn)=(H(nn)<=h3(nn)).*Ol_ex(nn)+...
          (vf_aw(nn)>wsp_const).*(H(nn)<=h4(nn)).*(H(nn)>h3(nn)).*spR_dgr(nn).*(1-vf_aw(nn))+...
          (vf_aw(nn)<=wsp_const).*(H(nn)<=h4(nn)).*(H(nn)>h3(nn)).*(spR_dgr(nn).*(1-vf_aw(nn)-Ol_ex(nn))+Ol_ex(nn)) +...             
          (1-vf_aw(nn)).*(H(nn)>h4(nn))    ;

v_ol(nW)=0;
Ol_ex(nW)=0;
v_sp=1-phi-v_ol;

vf_aw2=v_sp*wsp_const+phi;
DDWW=vf_aw2-vf_aw;
vf_aw=vf_aw2;

Ol_sys=v_ol-Ol_ex;   %system means the Olivine+free water=serpentine system  
wfC=Ol_sys*(wsp_const)/(1-wsp_const);  %volume fraction of potential combined water

w_ex=max(phi-wfC,0);   %extra free water beyond the serpentine system
V_sys=1-Ol_ex-w_ex;   %system means the Olivine+free water=serpentine system 

Cp=(phi.*Cpi.*(1-w1)*rho_i+phi.*Cpw.*w1*rho_w+v_sp*CpSp*rho_sp+v_ol*CpOl*rho_ol)/Ho;   %Cp after the thermal step

n1=(phi>0);
nxn=(phi==0);

w1(n1)=(vf_aw(n1)<=wsp_const).* (H(n1)>=h3(n1)) + ...
    (vf_aw(n1)>wsp_const).*((H(n1)-h1(n1))./(h2(n1)-h1(n1)).*((H(n1)<h2(n1)).*(H(n1)>h1(n1)))+(H(n1)>=h2(n1)));

w1(nxn)=1; 

 T_dL=(vf_aw<=wsp_const).* ((H<=h3).*(th0+H./Cp) + T_sp.*(H>h3).*(H<=h4)...
     +(T_sp+(H-h4)./Cp).*(H>h4));    

 T_dR=(vf_aw>wsp_const).* ((H<=h1).*(th0+H./Cp) + Tas.*(H>h1).*(H<=h2) + (Tas+(H-h2)./Cp).*(H>h2).*(H<=h3)...
     + T_sp.*(H>h3).*(H<=h4)  + (T_sp+(H-h4)./Cp).*(H>h4)); 

 T_d=T_dL+T_dR;

T=(T_d-Tus)/Tc;  

if i+1==LocbX
j=j+1;

    if isempty(zN)
        zN=1;
    
    else
        p_c(1:zN,j)=pfull(1:zN);
        pf_c(1:zN,j)=pfull(1:zN)+ps(1:zN);
        v_c(1:zN+1,j)=v_sfaceh(1:zN+1);
        qr_c(1:zN+1,j)=qrh(1:zN+1);
        stainR_c(:,j)=stainR;
    end

deVN(j)=deV;
deWN(j)=deW;
deW1N(j)=deW1;
DWN(:,j)=DDWW;
NNzone(j)=N_zone;
phi_c(:,j)=phi;
C_c(:,j)=C;

vf_sp_c(:,j)=vf_dgh;  %hydration degree, =1 when there's cold olivine but no free water
v_sp_c(:,j)=v_sp;
v_ol_c(:,j)=v_ol;

V_sys_c(:,j)=V_sys;
vf_aw_c(:,j)=vf_aw;
Ol_ex_c(:,j)=Ol_ex;
w_ex_c(:,j)=w_ex;


H_c(:,j)=H;
T_c(:,j)=T;
w_c(:,j)=w1;
ps_c(:,j)=psD;

nzN(j)=zN;

if N_zone>=1
    nzN2(j)=i_right(N_zone);
else
    nzN2(j)=1;
end
S_Phid=sum(delta^3*V_ra'*phi/Vo,1)-S_Phi0;


end

end

%%
set(groot,'defaultAxesTickLabelInterpreter','latex');  

%%

 rho_n1=rho_w*phi_c+rho_sp*(1-phi_c-v_ol_c)+rho_ol*v_ol_c;    %only in the current case: rho_w=rho_i
 m_t=Grid.V'*rho_n1*delta^3/M;

 figure(1)  %total mass conservation check
 plot(t2(1:j)/1e6,m_t(1:j))
 xlabel('Time since CAI (Myr)')

ylabel('M/M0')
%%

ntp=178; %determined manually in each case

V_wcore=Grid.V(1:ntp)'*phi_c(1:ntp,:)/sum(Grid.V(1:ntp));
V_olcore=Grid.V(1:ntp)'*v_ol_c(1:ntp,:)/sum(Grid.V(1:ntp));
V_serpcore=Grid.V(1:ntp)'*vf_sp_c(1:ntp,:)/sum(Grid.V(1:ntp));
V_fwcore=Grid.V(1:ntp)'*(phi_c(1:ntp,:).*w_c(1:ntp,:))/sum(Grid.V(1:ntp));

CoreTopR=Grid.xf(ntp+1)*delta;   %m

%% Conservation check for total water volume fraction
totW=Grid.V'*vf_aw_c/Vc;   %volume fraction of free (ice+liquid) + combined water
totSW=Grid.V'*DWN/Vc;
figure(5)
 plot(t2/1e6,totW)
 xlabel('Time since CAI (Myr)')

ylabel('vf_aw')

%% moment of inertia
Isphc=8/3*pi*(grxDm.^4*rho_n1*dx0)*delta;
IoMRc=Isphc./(m_t*M*R^2);


%%
figure(4)

plot(t2/1e6,IoMRc)
ylabel('Moment of inertia')
xlabel('Time since CAI (Myr)')

%%   figure_16
Z1=Tc*T_c+Tus;

figW = 600; %scw/6;
gapL = 55;  %gap on left boundary
gapR = 40;  %gap on right boundary
gapT = 30;  %gap on top boundary
gapB = 60;  %gap on bottom boundary
gapH=34;

gapW=70;

subW = (figW-gapL-gapR);
subH = 0.25*subW;
figH = 5*subH+gapT+gapB+4*gapH;
subx_a =  gapL/figW;
suby_a = (gapB+4*(gapH+subH))/figH;
subx_b =  gapL/figW;
suby_b = (gapB+3*(gapH+subH))/figH;
subx_c =  subx_a;
suby_c = (gapB+2*(gapH+subH))/figH;
subx_d =  subx_a;
suby_d = (gapB+1*(gapH+subH))/figH;
subx_e =  subx_a;
suby_e = (gapB)/figH;


f16=figure(16);
set(f16,'Position', [0, 0, figW, figH])

AA=t2/1e6;
nt0=size(AA,2);
nt=size(AA,2);

Za=Z1;
h16_1=subplot('position',[subx_a suby_a subW/figW subH/figH]);

bottom = min(min(Za));
top  = max(max(Za));

hold on
X=repmat(AA,Nx1,1);
Y=repmat(delta*grx0/1000,size(AA,2),1)';
levels=linspace(bottom,top,50);


contourf(X,Y,Za,levels,'edgecolor','none')
contour(X,Y,Z1,[Tas Tas],'r') 
contour(X,Y,Z1,[T_sp T_sp],'k') 
ax = gca;
ax.Layer = 'top';
c = colorbar;
ylabel('R (km)')
c.Label.String = 'T (K)';
c.Label.FontSize=12;
hold off
xticks([tf/1e6 1000 2000 3000 4000 5000])
yticks([delta*grx0(1)/1000 100 200 300 400 delta*grx0(end)/1000])
yticklabels({'0', '100' ,'200', '300', '400' ,'470'})
ax.FontSize = 11;
text(0,500,'(a)')


Zb=phi_c;
h16_2=subplot('position',[subx_b suby_b subW/figW subH/figH]);
hold on
X=repmat(AA,Nx1,1);
Y=repmat(delta*grx0/1000,size(AA,2),1)';
levels=linspace(0,1,50);
contourf(X,Y,Zb,levels,'edgecolor','none')
contour(X,Y,Z1,[Tas Tas],'r') 
contour(X,Y,Z1,[T_sp T_sp],'k') 
c = colorbar;
clim manual
clim([0 1]);
ylabel('R (km)')
c.Label.String = 'free water';
c.Label.FontSize=12;
hold off
text(0,500,'(b)')
xticks([tf/1e6 1000 2000 3000 4000 5000])
yticks([delta*grx0(1)/1000 100 200 300 400 delta*grx0(end)/1000])
yticklabels({'0', '100' ,'200', '300', '400' ,'470'})
ax = gca;
ax.FontSize = 11;
ax.Layer = 'top';


h16_3=subplot('position',[subx_c suby_c subW/figW subH/figH]);
Zc=v_ol_c;
hold on
X=repmat(AA,Nx1,1);
Y=repmat(delta*grx0/1000,size(AA,2),1)';
levels=linspace(0,1,50);
contourf(X,Y,Zc,levels,'edgecolor','none')
contour(X,Y,Z1,[T_sp T_sp],'k') 
c = colorbar;
clim auto
clim([0 1]);
ylabel('R (km)')
c.Label.String = 'Olivine';
c.Label.FontSize=12;
hold off
text(0,500,'(c)')
xticks([tf/1e6 1000 2000 3000 4000 5000])
yticks([delta*grx0(1)/1000 100 200 300 400 delta*grx0(end)/1000])
yticklabels({'0', '100' ,'200', '300', '400' ,'470'})
ax = gca;
ax.FontSize = 11; 
ax.Layer = 'top';


Zd=(1-phi_c-v_ol_c);
Zd(Zd<0)=0;
h16_4=subplot('position',[subx_d suby_d subW/figW subH/figH]);
hold on
X=repmat(AA,Nx1,1);
Y=repmat(delta*grx0/1000,size(AA,2),1)';
levels=linspace(0,1,50);
contourf(X,Y,Zd,levels,'edgecolor','none')
contour(X,Y,Z1,[T_sp T_sp],'k') 
ylabel('R (km)')
c = colorbar;
clim manual
clim([0 1]);
c.Label.String = 'Serpentine';
c.Label.FontSize=12;
hold off
xticks([tf/1e6 1000 2000 3000 4000 5000])
yticks([delta*grx0(1)/1000 100 200 300 400 delta*grx0(end)/1000])
yticklabels({'0', '100' ,'200', '300', '400' ,'470'})
text(0,500,'(d)')
ax = gca;
ax.FontSize = 11; 
ax.Layer = 'top';

h16_5=subplot('position',[subx_e suby_e 0.851*subW/figW subH/figH]);
hold on
plot(t2/1e6,1e2*(V_wcore+V_olcore),'LineWidth',1,'Color',[0.8500 0.3250 0.0980])
plot(t2/1e6,1e2*V_wcore,'LineWidth',1,'Color',[0 0.4470 0.7410]);
hold off


xticks([tf/1e6 1000 2000 3000 4000 5000])        
xlabel('Time since CAI (Myr)')
ylabel('Core Volume (%)')
ax = gca;
ax.FontSize = 11; 

axis([t2(1)/1e6 t2(end)/1e6 5 20])
text(800,15,'Olivine','FontSize',12)
text(700,9,'Free Water','FontSize',12)
text(3000,16.5,'Serpentine','FontSize',12)
text(0,21,'(e)')


%%  Movie  1
Myr2s=1e6*yr2s;
figW=1200;
figH=500;
linw=1;
nzN(nzN==0)=1;
V = VideoWriter(['DCZ',num2str(te/1e9),'_',num2str(R/1e3)]);

V.FrameRate = 30; 
    open(V); 
    nsp=6;

 for ii=1:j
 f12=figure(12);
 set(f12,'Position', [500, 0, figW, figH])
xx=ii;

 subplot(1,nsp,1)
   plot(Myr2s*qrc*qr_c(1:nzN(xx)+1,xx),delta/1e3*grxf(1:nzN(xx)+1),'b'...
       ,Myr2s*qrc*v_c(1:nzN(xx)+1,xx),delta/1e3*grxf(1:nzN(xx)+1),'r','LineWidth',linw)
xlabel('v_s r & q_r b (m/Myr)')
  ylabel('Radius (km)')
ylim([0 R/1e3])
xlim([-32 32])
text(-20,1.05*R/1e3,['t=',num2str(t2(xx)/1e6),' Myr'])


 subplot(1,nsp,2)
 plot(pc/1e6*p_c(1:nzN(xx),xx),delta/1e3*grx0(1:nzN(xx)),zeros(nzN(xx),1)...
     ,delta/1e3*grx0(1:nzN(xx)),'k--','LineWidth',linw)
 xlabel('p (MPa)')
ylim([0 R/1e3])
xlim([-5 5])


  subplot(1,nsp,3)
plot(T_SI,delta/1e3*grx0,'k--',T_sp*ones(Nx1,1),delta/1e3*grx0,'r--',Tc*T_c(:,xx)+Tus,delta/1e3*grx0,'b','LineWidth',linw)
 xlabel('T (K)')
ylim([0 R/1e3])
xlim([155 900])
ylabel('Radius (km)')

 subplot(1,nsp,4)
 plot(phi_c(:,xx),delta*grx0/1e3,phi_c(nzN(xx),xx),delta*grx0(nzN(xx))/1e3...
     ,'.','MarkerSize',10,'LineWidth',linw)
 xlabel('\phi')
axis([0 1 0 R/1e3])


 subplot(1,nsp,5)
  plot(vf_sp_c(:,xx),delta/1e3*grx0,'b','LineWidth',linw)
 xlabel('''degree of hydration''')
ylim([0 R/1e3])
xlim([0 1])


subplot(1,nsp,6)
  plot(v_ol_c(:,xx),delta/1e3*grx0,'b','LineWidth',linw)
 xlabel('olivine volume fraction')
ylim([0 R/1e3])
xlim([0 1])

   frame=getframe(gcf);
 writeVideo(V,frame);
 end

 close(V);

 %% write result
% addsp1='D:\Projects\Ceres compaction\compaction\data26\';
% addsp2='D:\Projects\Ceres compaction\compaction\data26\MOI_tY\';
% name1='cm1e-12_mu3.xlsx';  %_mu-2
% writematrix(qr_c,[addsp1,name1],'Sheet', 'qr_c', 'Range','A1');
% writematrix(v_c',[addsp1,name1],'Sheet', 'v_c', 'Range','A1');
% writematrix(p_c',[addsp1,name1],'Sheet', 'p_c', 'Range','A1');
% writematrix(phi_c',[addsp1,name1],'Sheet', 'phi_c', 'Range','A1');
% writematrix(vf_sp_c',[addsp1,name1],'Sheet', 'vf_sp_c', 'Range','A1');
% writematrix(v_ol_c',[addsp1,name1],'Sheet', 'v_ol_c', 'Range','A1');
% writematrix(T_c',[addsp1,name1],'Sheet', 'T_c', 'Range','A1');
% writematrix(H_c',[addsp1,name1],'Sheet', 'H_c', 'Range','A1');
% writematrix(C_c',[addsp1,name1],'Sheet', 'C_c', 'Range','A1');
% writematrix(w_c',[addsp1,name1],'Sheet', 'w_c', 'Range','A1');
% 
% 
% writematrix(t2',[addsp2,name1],'Sheet', 't', 'Range','A1');
% writematrix(IoMRc',[addsp2,name1],'Sheet', 'IC', 'Range','A1');
% writematrix(1e2*(V_wcore+V_olcore)',[addsp2,name1],'Sheet', 'olphi', 'Range','A1');
% writematrix(1e2*(V_wcore)',[addsp2,name1],'Sheet', 'phi', 'Range','A1');
% writematrix(1e2*(V_fwcore)',[addsp2,name1],'Sheet', 'fw', 'Range','A1');
