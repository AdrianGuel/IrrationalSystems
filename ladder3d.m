clear all; close all; clc;

w=0.001:1:130000;
C=0.1e-6;
L=1e-3;
w0=1/sqrt(L*C);
z0=sqrt(L/C);
mu = 0.7;

kd0=0:.1:10;
kp0=-z0*ones(1,length(kd0));
z=ones(1,length(w)); z2=ones(1,length(kd0));


fighdl = figure(); view(3); hold on;
%view([0,0,1]);
set(fighdl,'defaulttextInterpreter','latex') 
set(gcf,'color','w');
axis square
axis([-105 -80 -20 10 0 1])
grid on

xlabel('$k_p(\omega)$','FontSize', 16)
ylabel('$k_\eta(\omega)$','FontSize', 16)
zlabel('$\alpha$','FontSize', 16)

for mu=0.01:.1:1
[ d, ix ] = min( abs( kpx(w,mu,w0,z0))-80);
aux1=kdx(w,mu,w0,z0);
plot3(kpx(w,mu,w0,z0),aux1,z.*mu,'--b',kp0,kd0,z2.*mu,'r')%,[-100 -80],[5 5],[mu mu],'--b',[-80 -80],[8 aux1(ix)],[mu mu],'--b')
%plot(kpx(w,mu),kdx(w,mu),'b','LineWidth',3)
end

kpi=-105:0.1:-80;
kdi=-(z0/w0)*ones(1,length(kpi));
plot3(kpx(w,1,w0,z0),kdx(w,1,w0,z0),z.*1,'b','LineWidth',3)
plot3([-100 -100],[10 -20],[1 1],'r','LineWidth',3)
plot3(kpi,kdi,ones(1,length(kdi)).*1,'y','LineWidth',3)
plot3([-105 -100],[-20 -20],[1 1],'k--','LineWidth',3)
plot3([-105 -105],[-20 0],[1 1],'k--','LineWidth',3)
% kd0=0:.1:1;
% kp0=-ones(1,length(kd0));
% aux2=kpx(w,0.5);
% kpi=aux2(end):0.1:1;
% kdi=kpi;
% kdi(:)=-1/sqrt(3);
% plot3(kp0,kd0,ones(1,length(kd0)).*0.5,'r','LineWidth',3)
% plot3(kpi,kdi,ones(1,length(kdi)).*0.5,'y','LineWidth',3)



 kp=[-90 -99.5 -100.3];
 kd=[3 -0.007 0];
 syms q
for mu=0.01:.02:1
dl=vpasolve(fragility(kp(1),kd(1),q,mu,w0,z0)==0,q,[0 130000]);
if ~isempty(dl)
    d_02=kp(1)+z0;
    %d_inf2=kd(1)+1/sqrt(3);
    d_l2=sqrt((kpx(dl,mu,w0,z0)-kp(1))^2+(kdx(dl,mu,w0,z0)-kd(1))^2);
    radius = min(abs([d_02 d_l2])); 
    xCenter = kp(1);
    yCenter = kd(1);    
    theta = 0 : 0.01 : 2*pi;
    x = radius * cos(theta) + xCenter;
    y = radius * sin(theta) + yCenter;
    plot3(x,y,mu*ones(1,length(x)),'k--');
 %   arrow([xCenter yCenter],[x(end) y(end)],'Color','r','Length',9)
 %   text((xCenter+x(end))/2-0.1,(yCenter+y(end))/2-0.2,'$d_2$','Interpreter','Latex','Color','k','FontSize',10)
else
    d_02=kp(1)+z0;
    %d_inf2=kd(2)+1/sqrt(3);
    radius = min(abs([d_02])); 
    xCenter = kp(1);
    yCenter = kd(1);    
    theta = 0 : 0.01 : 2*pi;
    x = radius * cos(theta) + xCenter;
    y = radius * sin(theta) + yCenter;
    plot3(x,y,mu*ones(1,length(x)),'k--');
end
end
arrow3([kp(1) kd(1) 1],[kp(1) kd(1) .01])
plot3(kp(1),kd(1),1,'kp','LineWidth',2)
text(kp(1)+.01,kd(1),1,'$\mathbf{k}_4$','Interpreter','Latex','Color','k','FontSize',14)

function kpaux=kpx(w,mu,w0,z0)
kpaux=(1/2).*w0.^(-1).*z0.*csc(mu.*atan2(w,0)).*(w.*cos(mu.*atan2(w,0))+(-1).* ...
  w0.*(w0.^(-4).*(w.^2+(-4).*w0.^2).^2).^(1/4).*sin(mu.*atan2(w,0)+(-1/2) ...
  .*atan2(0,4+(-1).*w.^2.*w0.^(-2))));
end

function kpaux=kdx(w,mu,w0,z0)
kpaux=(-1/2).*(w.^2).^((-1/2).*mu).*w0.^(-1).*z0.*csc(mu.*atan2(w,0)).*(w+w0.* ...
  (w0.^(-4).*(w.^2+(-4).*w0.^2).^2).^(1/4).*sin((1/2).*atan2(0,4+(-1).* ...
  w.^2.*w0.^(-2))));
end


function Omegal=fragility(kp,kd,w,mu,w0,z0)
Omegal=(1/2).*w0.^(-4).*z0.*csc(mu.*atan2(w,0)).*((1/2).*w0.^2.*(cos(mu.*atan2( ...
  w,0))+w.*w0.*(w0.^(-4).*(w.^2+(-4).*w0.^2).^2).^(1/4).*((-1).*w.^2+4.* ...
  w0.^2).^(-1).*sin(mu.*atan2(w,0)+(-1/2).*atan2(0,4+(-1).*w.^2.*w0.^(-2)) ...
  )).*(w.*z0.*cot(mu.*atan2(w,0))+(-1).*w0.*(2.*kp+(w0.^(-4).*(w.^2+(-4).* ...
  w0.^2).^2).^(1/4).*z0.*csc(mu.*atan2(w,0)).*sin(mu.*atan2(w,0)+(-1/2).* ...
  atan2(0,4+(-1).*w.^2.*w0.^(-2)))))+(w.^2).^((-1/2).*mu).*(((-1)+mu).* ...
  w0.^3+w.^(-1).*(w.^2+(-4).*w0.^2).*(w0.^(-4).*(w.^2+(-4).*w0.^2).^2).^( ...
  -3/4).*(((-1)+mu).*w.^2+(-4).*mu.*w0.^2).*sin((1/2).*atan2(0,4+(-1).* ...
  w.^2.*w0.^(-2)))).*((-1).*kd+(-1/2).*(w.^2).^((-1/2).*mu).*w0.^(-1).* ...
  z0.*csc(mu.*atan2(w,0)).*(w+w0.*(w0.^(-4).*(w.^2+(-4).*w0.^2).^2).^(1/4) ...
  .*sin((1/2).*atan2(0,4+(-1).*w.^2.*w0.^(-2))))));
end