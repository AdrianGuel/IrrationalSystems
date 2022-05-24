clear all
close all
w=0.001:1:3900;
C=0.1e-6;
L=1e-3;
w0=1/sqrt(L*C);
z0=sqrt(L/C);
mu = 0.7;

kd0=-2:.1:2;
kp0=-z0*ones(1,length(kd0));

fighdl = figure(2);
%view([0,0,1]);
set(fighdl,'defaulttextInterpreter','latex') 
set(gcf,'color','w');

%axis([-1.5 1.5 -1.5 0.5])

subplot(2,2,[1,3])
 kpm=-100:1:-99;
 kdm=0.01*ones(length(kpm));
 H=area(kpm,kdm,-2);
 set(H(1),'FaceColor',[0.85 0.85 0.85],'LineStyle','none','ShowBaseLine','off');
 hold on
 H=area(kpx(w,mu,w0,z0),kdx(w,mu,w0,z0),-2);
 set(H(1),'FaceColor',[1 1 1],'LineStyle','none','ShowBaseLine','off');
 hold on
plot(kpx(w,mu,w0,z0),kdx(w,mu,w0,z0),'b','LineWidth',2)
hold on
plot(kp0,kd0,'r','LineWidth',2);

x = [.39 0.35];
y = [.4 0.4];
annotation('textarrow',x,y,'String','$CRB$','Interpreter','Latex','FontSize', 12)
x = [.27 0.24];
y = [.7 0.7];
annotation('textarrow',x,y,'String','$RRB$','Interpreter','Latex','FontSize', 12)


%pickup 4 points
%[kp,kd] = ginput(5);
 kp=[-99.8 -99.5 -100.3];
 kd=[0 -0.007 0];
hold on
plot(kp(1),kd(1),'kp','LineWidth',0.5)
text(kp(1)+.1,kd(1),'$\mathbf{k}_1$','Interpreter','Latex','Color','k','FontSize',14)
hold on
plot(kp(2),kd(2),'kp','LineWidth',1)
text(kp(2)+.1,kd(2),'$\mathbf{k}_2$','Interpreter','Latex','Color','k','FontSize',14)
hold on
plot(kp(3),kd(3),'kp','LineWidth',1)
text(kp(3)+.1,kd(3),'$\mathbf{k}_3$','Interpreter','Latex','Color','k','FontSize',14)
hold on

syms x
dl=vpasolve(fragility(kp(1),kd(1),x,mu,w0,z0)==0,x,[0 900]);
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
    plot(x,y,'k--');
    hold on
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
    plot(x,y,'k--')
    hold on
end

txt = 'Stability region';
text(-99.5,0,txt);
xlabel('$k_p(\omega)$','FontSize', 16)
ylabel('$k_\eta(\omega)$','FontSize', 16)
axis([-100.5 -99 -0.01 0.01])
axis square

axes('Position',[.35 .6 .15 .15])
box on
plot(x,y,'k--')
hold on
plot(kp(1),kd(1),'kp','LineWidth',0.5)
hold on
w=740:.1:820;
plot(kpx(w,mu,w0,z0),kdx(w,mu,w0,z0),'b','LineWidth',2)

w=0.001:1:3900;
subplot(2,2,2)
plot(kpx(w,mu,w0,z0),real(DSkp(w,mu,w0,z0,kdx(w,mu,w0,z0))),'k--','LineWidth',1)
axis('square')
xlabel('$k_p(\omega)$','Interpreter','Latex','FontSize', 16)
ylabel('$\mathcal{S}_p$','Interpreter','Latex','FontSize', 16)
axis([-100 -99 -2800 -2600]);

subplot(2,2,4)
plot(kdx(w,mu,w0,z0),real(DSkd(w,mu,w0,z0,kdx(w,mu,w0,z0))),'k--','LineWidth',1)
axis('square')
xlabel('$k_\eta(\omega)$','Interpreter','Latex','FontSize', 16)
ylabel('$\mathcal{S}_\eta$','Interpreter','Latex','FontSize', 16)
axis([-0.007 0 -15e5 0]);



%[kp,kd] = ginput(3);
% %Simulation impulse response
t = 0:.00001:.1;
result1=euler_inversion(@(s) (2.*(kp(1)+kd(1).*s.^mu).*w0.*(2.*kp(1).*w0+2.*kd(1).*s.^mu.*w0+(s+(4+s.^2.*w0.^(-2)) ...
  .^(1/2).*w0).*z0).^(-1))*(1/s), t);
result2=euler_inversion(@(s) (2.*(kp(2)+kd(2).*s.^mu).*w0.*(2.*kp(2).*w0+2.*kd(2).*s.^mu.*w0+(s+(4+s.^2.*w0.^(-2)) ...
  .^(1/2).*w0).*z0).^(-1))*(1/s), t);
result3=euler_inversion(@(s) (2.*(kp(3)+kd(3).*s.^mu).*w0.*(2.*kp(3).*w0+2.*kd(3).*s.^mu.*w0+(s+(4+s.^2.*w0.^(-2)) ...
  .^(1/2).*w0).*z0).^(-1))*(1/s), t);
%result4=euler_inversion(@(s) ((kp(4)+kd(4)*s^mu)/(sqrt(s+1)+(kp(4)+kd(4)*s^mu))), t);

figure
subplot(3,1,1)
plot(t,result1,'b','LineWidth',2)
axis([0 0.1 -600 0]);
%axis('square')
%xlabel('Time (sec.)','Interpreter','Latex','FontSize', 14)
ylabel('$y(t)$','Interpreter','Latex','FontSize', 16)
str = sprintf('$k_1=(k_p,k_\\eta)=$ (%.2f,%.2f)',kp(1),kd(1));
title(str,'Interpreter','Latex')

subplot(3,1,2)
plot(t,result2,'b','LineWidth',2)
axis([0 0.013 -1e13 1e13]);
%axis('square')
%xlabel('Time (sec.)','Interpreter','Latex','FontSize', 14)
ylabel('$y(t)$','Interpreter','Latex','FontSize', 16)
str = sprintf('$k_2=(k_p,k_\\eta)=$ (%.2f,%.2f)',kp(2),kd(2));
title(str,'Interpreter','Latex')
subplot(3,1,3)
plot(t,result3,'b','LineWidth',2)
axis([0 0.0409 -11e13 1e3]);
%axis('square')
xlabel('Time (sec.)','Interpreter','Latex','FontSize', 14)
ylabel('$y(t)$','Interpreter','Latex','FontSize', 16)
str = sprintf('$k_3=(k_p,k_\\eta)=$ (%.2f,%.2f)',kp(3),kd(3));
title(str,'Interpreter','Latex')
%subplot(2,2,4)
%plot(t,result4)
set(gcf,'color','w');

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

function Dkpaux=DSkp(w,mu,w0,z0,kd)
Dkpaux=(-2).*w0.*(2.*kd.*mu.*(sqrt(-1).*w).^((-1)+mu).*w0+(1+sqrt(-1).*w.*(4+( ...
  -1).*w.^2.*w0.^(-2)).^(-1/2).*w0.^(-1)).*z0).^(-1);
end

function Dkdaux=DSkd(w,mu,w0,z0,kd)
Dkdaux=(-2).*(sqrt(-1).*w).^mu.*w0.*(2.*kd.*mu.*(sqrt(-1).*w).^((-1)+mu).*w0+( ...
  1+sqrt(-1).*w.*(4+(-1).*w.^2.*w0.^(-2)).^(-1/2).*w0.^(-1)).*z0).^(-1);
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