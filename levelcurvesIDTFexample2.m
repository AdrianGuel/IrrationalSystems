clear all; close all; clc;


fighdl = figure(); view(3); hold on;
%view([0,0,1]);
set(fighdl,'defaulttextInterpreter','latex') 
set(gcf,'color','w');
axis square
grid on
xlabel('$k_p(\omega)$','FontSize', 16)
ylabel('$k_\eta(\omega)$','FontSize', 16)
zlabel('$\alpha$','FontSize', 16)
axis([-1.5 0 -1.5 0.5 0 0.55]); 
kd0=0:1:3; w=0.00001:.01:800;
kp0=-ones(1,length(kd0));

z=ones(1,length(w)); z2=ones(1,length(kd0));


for mu=0.01:.01:0.5
%[ d, ix ] = min( abs( kpx(w,mu))-0.2);
aux1=kdx(w,mu);
plot3(kpx(w,mu),aux1,z.*mu,'--b',kp0,kd0,z2.*mu,'r')%,[-1 0],[.5 .5],[mu mu],'--b',[0 0],[0.5 aux1(ix)],[mu mu],'--b')
%plot(kpx(w,mu),kdx(w,mu),'b','LineWidth',3)
end
%aux1=kdx(w,0.6);
%plot3(kpx(w,0.6),aux1,z.*0.6,'b','LineWidth',3)
%X = [-1.5 -1 -1 -1.5 ];
%Y = [-1.5 -1.5  -1/sqrt(3) -1/sqrt(3)];
%Z = [0.5 0.5 0.5 0.5];
%fill3(X,Y,Z,[0.95 0.95 0.95])
kpi=-5:0.1:-1;
kdi=kpi;
kdi(:)=-1/sqrt(3);
kd0=-2:.1:-1/sqrt(3);
kp0=-ones(1,length(kd0));
plot3(kpx(w,0.5),kdx(w,0.5),z.*0.5,'b','LineWidth',3)
plot3(kp0,kd0,ones(1,length(kd0)).*0.5,'r','LineWidth',3)
plot3(kpi,kdi,ones(1,length(kdi)).*0.5,'y','LineWidth',3)
plot3([-1 0],[0.5 0.5],[0.5 0.5],'b','LineWidth',3)
plot3([0 0],[0.5 -1/sqrt(3)],[0.5 0.5],'b','LineWidth',3)
kd0=0:.1:1;
kp0=-ones(1,length(kd0));
aux2=kpx(w,0.5);
kpi=aux2(end):0.1:1;
kdi=kpi;
kdi(:)=-1/sqrt(3);
plot3(kp0,kd0,ones(1,length(kd0)).*0.5,'r','LineWidth',3)
plot3(kpi,kdi,ones(1,length(kdi)).*0.5,'y','LineWidth',3)
plot3([-1.5 -1],[-1.5 -1.5],[0.5 0.5],'k--','LineWidth',3)
plot3([-1.5 -1.5],[-1.5 -1/sqrt(3)],[0.5 0.5],'k--','LineWidth',3)

kp=[-1.25 -0.80 -1.25 -0.68 -0.92];
kd=[-1.06 -0.13 0 -1.12 -0.44];

%plot3(kp(1),kd(1),0.5,'kp','LineWidth',2)
text(kp(2)-0.1,kd(2)+0.2,0.4,'$\mathbf{k}_2$','Interpreter','Latex','Color','k','FontSize',14)


% d_01=kp(1)+1;
% d_inf1=kd(1)+1/sqrt(3);
% radius = min(abs([d_01 d_inf1])); 
% xCenter = kp(1);
% yCenter = kd(1);
% theta = 0 : 0.01 : 2*pi;
% x = radius * cos(theta) + xCenter;
% y = radius * sin(theta) + yCenter;
% plot3(x,y,0.5*ones(1,length(x)),'k--');
%hold on
%arrow([xCenter yCenter],[x(end) y(end)],'Color','r','Length',9)
%text((xCenter+x(end))/2-0.1,(yCenter+y(end))/2-0.2,'$d_1$','Interpreter','Latex','Color','k','FontSize',10)
hold on
for mu=0.1:.1:0.4
syms x
dl=vpasolve(fragilityIDTFexample2(kp(2),kd(2),x,mu)==0,x,[0 2])
if ~isempty(dl)
    d_02=kp(2)+1
    d_inf2=kd(2)+1/sqrt(3);
    d_l2=sqrt((kpx(dl,mu)-kp(2))^2+(kdx(dl,mu)-kd(2))^2)
    radius = round(min(abs([d_02 d_inf2 d_l2]))*100)/100; 
    xCenter = kp(2);
    yCenter = kd(2);    
    theta = 0 : 0.01 : 2*pi;
    x = radius * cos(theta) + xCenter;
    y = radius * sin(theta) + yCenter;
    plot3(x,y,mu*ones(1,length(x)),'k--');
    plot3(x,y,mu*ones(1,length(x)),'k--')
 %   hold on
 %   arrow([xCenter yCenter],[x(end) y(end)],'Color','r','Length',9)
 %   text((xCenter+x(end))/2-0.1,(yCenter+y(end))/2-0.2,'$d_2$','Interpreter','Latex','Color','k','FontSize',10)
else
    d_02=kp(2)+1;
    d_inf2=kd(2)+1/sqrt(3);
    radius = min(abs([d_02 d_inf2])); 
    xCenter = kp(2);
    yCenter = kd(2);    
    theta = 0 : 0.01 : 2*pi;
    x = radius * cos(theta) + xCenter;
    y = radius * sin(theta) + yCenter;
    plot3(x,y,mu*ones(1,length(x)),'k--');
  %  hold on
  %  arrow([xCenter yCenter],[x(end) y(end)],'Color','r','Length',9)
  %  text((xCenter+x(end))/2-0.1,(yCenter+y(end))/2-0.2,'$d_2$','Interpreter','Latex','Color','k','FontSize',10)
    
end
end
arrow3([kp(2) kd(2) 0.4],[kp(2) kd(2) 0.1])
hold off


fighdl = figure(2); view(3); hold on;
set(fighdl,'defaulttextInterpreter','latex') 
set(gcf,'color','w');
axis square
xlabel('$k_p(\omega)$','FontSize', 16)
ylabel('$k_\eta(\omega)$','FontSize', 16)
zlabel('$\alpha$','FontSize', 16)
axis([-4 0 -1 2.5 -1 0.1]); 
kd0=-1:1:0; w=0.00001:.001:100;
kp0=zeros(1,length(kd0));
kpi=-5:0.1:0.8;
kdi=kpi;
kdi(:)=-1/sqrt(3);
z=ones(1,length(w)); z2=ones(1,length(kd0));
for mu=-1:.01:-0.01
[ d, ix ] = min(abs(kdx(w,mu)-2.5));
aux1=kdx(w,mu);
aux2=kpx(w,mu);
plot3(kpx(w,mu),aux1,z.*mu,'--b',kd0,kp0,z2.*mu,'r',[0 0],[0 2.5],[mu mu],'--b',[0 aux2(ix)],[2.5 aux1(ix)],[mu mu],'--b')
%plot(kpx(w,mu),kdx(w,mu),'b','LineWidth',3)
end

hold off


function kpaux=kpx(w,mu)
kpaux=(-1).*(1+9.*w.^2).^(-1/4).*csc(mu.*angle(sqrt(-1).*w)).*(cos((1/2).* ...
  angle(1+(sqrt(-1)*3).*w)).^2+sin((1/2).*angle(1+(sqrt(-1)*3).*w)).^2).^( ...
  -1).*((-1).*w.*cos((1/2).*angle(1+(sqrt(-1)*3).*w)).*cos(mu.*angle(sqrt( ...
  -1).*w))+(-1).*(1+4.*w.^2).^(1/4).*cos((1/2).*angle(1+(sqrt(-1)*3).*w)) ...
  .*cos(mu.*angle(sqrt(-1).*w)).*sin((1/2).*angle(1+(sqrt(-1)*2).*w))+(1+ ...
  4.*w.^2).^(1/4).*cos((1/2).*angle(1+(sqrt(-1)*2).*w)).*cos(mu.*angle( ...
  sqrt(-1).*w)).*sin((1/2).*angle(1+(sqrt(-1)*3).*w))+(1+4.*w.^2).^(1/4).* ...
  cos((1/2).*angle(1+(sqrt(-1)*2).*w)).*cos((1/2).*angle(1+(sqrt(-1)*3).* ...
  w)).*sin(mu.*angle(sqrt(-1).*w))+w.*sin((1/2).*angle(1+(sqrt(-1)*3).*w)) ...
  .*sin(mu.*angle(sqrt(-1).*w))+(1+4.*w.^2).^(1/4).*sin((1/2).*angle(1+( ...
  sqrt(-1)*2).*w)).*sin((1/2).*angle(1+(sqrt(-1)*3).*w)).*sin(mu.*angle( ...
  sqrt(-1).*w)));
end

function kdaux=kdx(w,mu)
kdaux=(-1).*(w.^2).^((-1/2).*mu).*(1+9.*w.^2).^(-1/4).*csc(mu.*angle(sqrt(-1) ...
  .*w)).*(w.*cos((1/2).*angle(1+(sqrt(-1)*3).*w))+(1+4.*w.^2).^(1/4).*cos( ...
  (1/2).*angle(1+(sqrt(-1)*3).*w)).*sin((1/2).*angle(1+(sqrt(-1)*2).*w))+( ...
  -1).*(1+4.*w.^2).^(1/4).*cos((1/2).*angle(1+(sqrt(-1)*2).*w)).*sin((1/2) ...
  .*angle(1+(sqrt(-1)*3).*w))).*(cos((1/2).*angle(1+(sqrt(-1)*3).*w)).^2+ ...
  sin((1/2).*angle(1+(sqrt(-1)*3).*w)).^2).^(-1);
end

function dskpaux=DSkp(w,mu,kp,kd)
dskpaux=(-1).*(1+(1+(sqrt(-1)*2).*w).^(-1/2)+(3/2).*(kp+kd.*(sqrt(-1).*w).^mu).* ...
  (1+(sqrt(-1)*3).*w).^(-1/2)+kd.*mu.*(1+(sqrt(-1)*3).*w).^(1/2).*(sqrt( ...
  -1).*w).^((-1)+mu)).^(-1).*(1+(sqrt(-1)*3).*w).^(1/2);
end

function dskdaux=DSkd(w,mu,kp,kd)
dskdaux=(-1).*(1+(1+(sqrt(-1)*2).*w).^(-1/2)+(3/2).*(kp+kd.*(sqrt(-1).*w).^mu).* ...
  (1+(sqrt(-1)*3).*w).^(-1/2)+kd.*mu.*(1+(sqrt(-1)*3).*w).^(1/2).*(sqrt( ...
  -1).*w).^((-1)+mu)).^(-1).*(1+(sqrt(-1)*3).*w).^(1/2).*(sqrt(-1).*w) ...
  .^mu;
end