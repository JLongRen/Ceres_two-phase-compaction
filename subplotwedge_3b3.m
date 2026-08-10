function subplotwedge_3b3(theta,Nx,Np,Field,Tcrt)

fontsz=15;
     xx=Field.phi(Field.b,:);
     water_1=Field.phi.*(Field.Tc>=Tcrt);
     fw1=water_1(Field.b,:);
     fol1=Field.v_ol(Field.b,:);

r=1;
x0=0;
y0=0;
linw=1;
   a1 = pi/2-theta/2;  
   a2 = a1 + theta;
   t = linspace(a1,a2);
   x = x0 + r*cos(t);
   y = y0 + r*sin(t);  %arch


   hold on
   
   plot([x0,x,x0],[y0,y,y0],'k','LineWidth',linw)
   pf1=fill([x0,x,x0],[y0,y,y0],[17, 119, 51]/255);   %col.green
   pf1.LineStyle = 'none';

   ND1=0:199;
   ND2=1:200;
   NDD1=repmat(ND1',1,Np);
   NDD2=repmat(ND2',1,Np);
   scale1=linspace(0,1,Np);
   scale2=linspace(1,0,Np);

   tt2=a1+theta*scale1'*xx;
   tt3=a1+theta*scale2'*xx;
    archB_x=x0+r*NDD1/200.*cos(tt2');  %vector for the bottom arch(x)
    archB_y=y0+r*NDD1/200.*sin(tt2');
    archT_x=x0+r*NDD2/200.*cos(tt3');
    archT_y=y0+r*NDD2/200.*sin(tt3');


    tw1=a1+theta*scale1'*fw1; 
    tw2=a1+theta*scale2'*fw1; 
    archBw_x=x0+r*NDD1/200.*cos(tw1');  %vector for the bottom arch(x)
    archBw_y=y0+r*NDD1/200.*sin(tw1');
    archTw_x=x0+r*NDD2/200.*cos(tw2');
    archTw_y=y0+r*NDD2/200.*sin(tw2');



    tol1=a1+theta*scale1'*fol1; 
    tol2=a1+theta*scale2'*fol1; 

    archBol_x=x0-r*NDD1/200.*cos(tol1');
    archBol_y=y0+r*NDD1/200.*sin(tol1');
    archTol_x=x0-r*NDD2/200.*cos(tol2');
    archTol_y=y0+r*NDD2/200.*sin(tol2');

for i=1:Nx

      pf4=fill([archBol_x(i,:),archTol_x(i,:)],[archBol_y(i,:),archTol_y(i,:)],[221, 204, 119]/255);
      pf4.LineStyle = 'none';

      pf3=fill([archB_x(i,:),archT_x(i,:)],[archB_y(i,:),archT_y(i,:)],[1, 1, 1]);
      pf3.LineStyle = 'none';

      pf2=fill([archBw_x(i,:),archTw_x(i,:)],[archBw_y(i,:),archTw_y(i,:)],[136, 204, 238]/255);%col.blue  [86, 180, 233]
      pf2.LineStyle = 'none';

end

%tick bars
t2=zeros(4,100);
x_b1=zeros(4,100);
y_b1=zeros(4,100);

t2r=zeros(4,100);
x_b1r=zeros(4,100);
y_b1r=zeros(4,100);

Rx=[100 200 300 400];
for i=1:4    
    t2(i,:)= linspace(a2,a2+1/30*(a2-a1)*400/Rx(i),100);
   x_b1(i,:) = x0 + Rx(i)/470*r*cos(t2(i,:));
   y_b1(i,:) = y0 + Rx(i)/470*r*sin(t2(i,:));  %arch

   t2r(i,:)= linspace(a1,a1-1/30*(a2-a1)*400/Rx(i),100);
   x_b1r(i,:) = x0 + Rx(i)/470*r*cos(t2r(i,:));
   y_b1r(i,:) = y0 + Rx(i)/470*r*sin(t2r(i,:));  %arch
end

for i=1:4
   plot(x_b1(i,:),y_b1(i,:),'k','LineWidth',linw)
   plot(x_b1r(i,:),y_b1r(i,:),'k','LineWidth',linw)
end

t3=zeros(5,100);
x_b2=zeros(4,100);
y_b2=zeros(4,100);
t3r=zeros(5,100);
x_b2r=zeros(4,100);
y_b2r=zeros(4,100);
Rx2=[50 150 250 350 450];
for i=1:5    
    t3(i,:)= linspace(a2,a2+1/40*(a2-a1)*400/Rx2(i),100);
   x_b2(i,:) = x0 + Rx2(i)/470*r*cos(t3(i,:));
   y_b2(i,:) = y0 + Rx2(i)/470*r*sin(t3(i,:));  %arch

       t3r(i,:)= linspace(a1,a1-1/40*(a2-a1)*400/Rx2(i),100);
   x_b2r(i,:) = x0 + Rx2(i)/470*r*cos(t3r(i,:));
   y_b2r(i,:) = y0 + Rx2(i)/470*r*sin(t3r(i,:));  %arch
end

for i=1:5
   plot(x_b2(i,:),y_b2(i,:),'k','LineWidth',linw)
   plot(x_b2r(i,:),y_b2r(i,:),'k','LineWidth',linw)
end


x_b4 =x_b1(1,:)-(x0 + Rx(1)/470*r*cos(a2));
y_b4 =y_b1(1,:)-(y0 + Rx(1)/470*r*sin(a2));
x_b4r =x_b1r(1,:)-(x0 + Rx(1)/470*r*cos(a1));
y_b4r =y_b1r(1,:)-(y0 + Rx(1)/470*r*sin(a1));
plot(x_b4,y_b4,'k','LineWidth',linw)
plot(x_b4r,y_b4r,'k','LineWidth',linw)

% xtext0=[-0.07 -0.22 -0.29 -0.38 -0.47 -0.515];
% ytext0=[-0.01 0.181 0.372 0.563 0.758 0.9];
% 
% xtext=linspace(1.3,1.15,6).*xtext0;
% ytext=0.95*ytext0;
xtext0=-0.06+[-0.07 -0.22 -0.29 -0.38 -0.47 -0.515];
ytext0=[-0.01 0.181 0.372 0.563 0.758 0.9];

xtext=linspace(1.3,1.15,6).*xtext0;
ytext=0.95*ytext0;
text(xtext(1),ytext(1),'0','interpreter','latex','FontSize',fontsz)
text(xtext(2),ytext(2),'100','interpreter','latex','FontSize',fontsz)
text(xtext(3),ytext(3),'200','interpreter','latex','FontSize',fontsz)
text(xtext(4),ytext(4),'300','interpreter','latex','FontSize',fontsz)
text(xtext(5),ytext(5),'400','interpreter','latex','FontSize',fontsz)
text(xtext(6),ytext(6),'470','interpreter','latex','FontSize',fontsz)
   hold off

   axis equal

set(gca,'xtick',[])
set(gca,'ytick',[])
box off
set(get(gca, 'XAxis'), 'Visible', 'off');
set(get(gca, 'YAxis'), 'Visible', 'off');