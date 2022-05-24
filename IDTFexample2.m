close all
clear all
mu=0.4;
w=0.00001:1:1000;

kd0=-5:1:3;
kp0=kd0;
kp0(:)=-1;

% kpi=-1.5:0.1:0.8;
% kdi=kpi;
% kdi(:)=-1/sqrt(3);

subplot(2,2,[1,3])
%kpm=-2:0.5:-1;
%kdm=zeros(length(kpm));
%kdm(:)=-1/sqrt(3);
%H=area(kpm,kdm,-2);
%set(H,'FaceColor',[0.85 0.85 0.85],'LineStyle','none','ShowBaseLine','off');
%hold on
kpm=-1:0.5:1;
kdm=zeros(length(kpm));
kdm(:)=0.5;
H2=area(kpm,kdm,-1.5);
set(H2,'FaceColor',[0.85 0.85 0.85],'LineStyle','none','ShowBaseLine','off');
hold on
w=0.00001:.01:50;
H3=area(kpx(w,mu),kdx(w,mu),-1.5);
set(H3,'FaceColor',[1 1 1],'LineStyle','none','ShowBaseLine','off');
hold on
w=0.00001:.01:100;
plot(kpx(w,mu),kdx(w,mu),'b','LineWidth',2)
hold on
plot(kp0,kd0,'r','LineWidth',2)
hold on
%plot(kpi,kdi,'y','LineWidth',2)
axis('square')
xlabel('$k_p(\omega)$','Interpreter','Latex','FontSize', 16)
ylabel('$k_\eta(\omega)$','Interpreter','Latex','FontSize', 16)
axis([-1.5 0 -1.5 0.5]);
txt = 'Stability Region';
text(-0.7,0.1,txt,'FontSize',8)
x = [.35 0.29];
y = [.53 0.53];
annotation('textarrow',x,y,'String','$CRB$','Interpreter','Latex','FontSize', 12)
x = [.3 0.25];
y = [.7 0.7];
annotation('textarrow',x,y,'String','$RRB$','Interpreter','Latex','FontSize', 12)
%x = [.4 0.4];
%y = [.45 0.5];
%annotation('textarrow',x,y,'String','$IRB$','Interpreter','Latex','FontSize', 12)

%pickup 4 points
%[kp,kd] = ginput(5);
kp=[-1.4 -0.80 -1.25 -0.68 -0.92];
kd=[-1.1 -0.13 0 -1.12 -0.44];
hold on
plot(kp(1),kd(1),'kp','LineWidth',2)
text(kp(1)-0.1,kd(1)+0.2,'$\mathbf{k}_1$','Interpreter','Latex','Color','k','FontSize',14)
hold on
plot(kp(2),kd(2),'kp','LineWidth',2)
text(kp(2)-0.1,kd(2)+0.2,'$\mathbf{k}_2$','Interpreter','Latex','Color','k','FontSize',14)
hold on
plot(kp(3),kd(3),'kp','LineWidth',2)
text(kp(3)-0.1,kd(3)+0.2,'$\mathbf{k}_3$','Interpreter','Latex','Color','k','FontSize',14)
hold on
plot(kp(4),kd(4),'kp','LineWidth',2)
text(kp(4)-0.1,kd(4)+0.2,'$\mathbf{k}_4$','Interpreter','Latex','Color','k','FontSize',14)
hold on
plot(kp(5),kd(5),'kp','LineWidth',2)
text(kp(5)-0.2,kd(5),'$\mathbf{k}_5$','Interpreter','Latex','Color','k','FontSize',14)
% hold on
% d_01=kp(1)+1;
% d_inf1=kd(1)+1/sqrt(3);
% radius = min(abs([d_01 d_inf1])); 
% xCenter = kp(1);
% yCenter = kd(1);
% theta = 0 : 0.01 : 2*pi;
% x = radius * cos(theta) + xCenter;
% y = radius * sin(theta) + yCenter;
% plot(x,y,'k--');
%hold on
%arrow([xCenter yCenter],[x(end) y(end)],'Color','r','Length',9)
%text((xCenter+x(end))/2-0.1,(yCenter+y(end))/2-0.2,'$d_1$','Interpreter','Latex','Color','k','FontSize',10)

hold on
syms x
dl=vpasolve(fragilityIDTFexample2(kp(2),kd(2),x,mu)==0,x,[0 2]);
if ~isempty(dl)
    d_02=kp(2)+1;
    d_inf2=kd(2)+1/sqrt(3);
    d_l2=sqrt((kpx(dl,mu)-kp(2))^2+(kdx(dl,mu)-kd(2))^2);
    radius = round(min(abs([d_02 d_inf2 d_l2]))*100)/100; 
    xCenter = kp(2);
    yCenter = kd(2);    
    theta = 0 : 0.01 : 2*pi;
    x = radius * cos(theta) + xCenter;
    y = radius * sin(theta) + yCenter;
    plot(x,y,'k--');
 %   hold on
 %   arrow([xCenter yCenter],[x(end) y(end)],'Color','r','Length',9)
 %   text((xCenter+x(end))/2-0.1,(yCenter+y(end))/2-0.2,'$d_2$','Interpreter','Latex','Color','k','FontSize',10)
%else
%    d_02=kp(2)+1;
%    d_inf2=kd(2)+1/sqrt(3);
%    radius = min(abs([d_02 d_inf2])); 
%    xCenter = kp(2);
%    yCenter = kd(2);    
%    theta = 0 : 0.01 : 2*pi;
%    x = radius * cos(theta) + xCenter;
%    y = radius * sin(theta) + yCenter;
%    plot(x,y,'k--');
  %  hold on
  %  arrow([xCenter yCenter],[x(end) y(end)],'Color','r','Length',9)
  %  text((xCenter+x(end))/2-0.1,(yCenter+y(end))/2-0.2,'$d_2$','Interpreter','Latex','Color','k','FontSize',10)
    
end

w=0.0001:0.01:20;
subplot(2,2,2)
plot(kpx(w,mu),real(DSkp(w,mu,kpx(w,mu),kdx(w,mu))),'k--')
axis('square')
xlabel('$k_p(\omega)$','Interpreter','Latex','FontSize', 16)
ylabel('$\mathcal{S}_p$','Interpreter','Latex','FontSize', 16)
axis([-1.5 0 -80 1]);
x = [.8 0.7];
y = [.41 0.37];
annotation('textarrow',x,y,'Color','red')


subplot(2,2,4)
plot(kdx(w,mu),real(DSkd(w,mu,kpx(w,mu),kdx(w,mu))),'k--')
axis('square')
xlabel('$k_\eta(\omega)$','Interpreter','Latex','FontSize', 16)
ylabel('$\mathcal{S}_\eta$','Interpreter','Latex','FontSize', 16)
axis([-1.5 0.5 -16 2]);
x = [.7 0.8];
y = [.91 0.89];
annotation('textarrow',x,y,'Color','red')
text(-0.3,-2,'$\omega$','Interpreter','Latex','Color','red','FontSize',14)
text(-0.3,26,'$\omega$','Interpreter','Latex','Color','red','FontSize',14)
set(gcf,'color','w');

%Simulation impulse response
t = 0:.1:100;
result1=euler_inversion(@(s) ((sqrt(3*s+1))*(kp(1)+kd(1)*s^mu)/(s+sqrt(2*s+1)+(kp(1)+kd(1)*s^mu)*(sqrt(3*s+1))))*(1/s), t);
result2=euler_inversion(@(s) ((sqrt(3*s+1))*(kp(2)+kd(2)*s^mu)/(s+sqrt(2*s+1)+(kp(2)+kd(2)*s^mu)*(sqrt(3*s+1))))*(1/s), t);
result3=euler_inversion(@(s) ((sqrt(3*s+1))*(kp(3)+kd(3)*s^mu)/(s+sqrt(2*s+1)+(kp(3)+kd(3)*s^mu)*(sqrt(3*s+1))))*(1/s), t);
result4=euler_inversion(@(s) ((sqrt(3*s+1))*(kp(4)+kd(4)*s^mu)/(s+sqrt(2*s+1)+(kp(4)+kd(4)*s^mu)*(sqrt(3*s+1))))*(1/s), t);
result5=euler_inversion(@(s) ((sqrt(3*s+1))*(kp(5)+kd(5)*s^mu)/(s+sqrt(2*s+1)+(kp(5)+kd(5)*s^mu)*(sqrt(3*s+1))))*(1/s), t);

figure
subplot(5,1,1)
plot(t,result1,'b','LineWidth',2)
%axis([0 50 0 4]);
%axis('square')
%xlabel('Time (sec.)','Interpreter','Latex','FontSize', 16)
ylabel('$y(t)$','Interpreter','Latex','FontSize', 16)
str = sprintf('${k}_1=(k_p,k_\\eta)=$ (%.2f,%.2f)',kp(1),kd(1));
title(str,'Interpreter','Latex')

subplot(5,1,2)
plot(t,result2,'b','LineWidth',2)
%axis('square')
%xlabel('Time (sec.)','Interpreter','Latex','FontSize', 16)
ylabel('$y(t)$','Interpreter','Latex','FontSize', 16)
str = sprintf('${k}_2=(k_p,k_\\eta)=$ (%.2f,%.2f)',kp(2),kd(2));
title(str,'Interpreter','Latex')
%axis([0 35 -4.5 0]);

subplot(5,1,3)
plot(t,result3,'b','LineWidth',2)
%axis('square')
%xlabel('Time (sec.)','Interpreter','Latex','FontSize', 16)
ylabel('$y(t)$','Interpreter','Latex','FontSize', 16)
str = sprintf('${k}_3=(k_p,k_\\eta)=$ (%.2f,%.2f)',kp(3),kd(3));
title(str,'Interpreter','Latex')
xlim([0 35])
%axis([0 35 -2e10 0]);
subplot(5,1,4)
plot(t,result4,'b','LineWidth',2)
%axis('square')
%xlabel('Time (sec.)','Interpreter','Latex','FontSize', 16)
ylabel('$y(t)$','Interpreter','Latex','FontSize', 16)
str = sprintf('${k}_4=(k_p,k_\\eta)=$ (%.2f,%.2f)',kp(4),kd(4));
title(str,'Interpreter','Latex')
%axis([0 35 0 2000]);
subplot(5,1,5)
plot(t,result4,'b','LineWidth',2)
%axis('square')
%axis([0 35 0 2000]);
xlabel('Time (sec.)','Interpreter','Latex','FontSize', 16)
ylabel('$y(t)$','Interpreter','Latex','FontSize', 16)
str = sprintf('${k}_5=(k_p,k_\\eta)=$ (%.2f,%.2f)',kp(5),kd(5));
title(str,'Interpreter','Latex')
set(gcf,'color','w');
% 

%figure
%w=0.00001:.001:2;
%plot(kpx(w,mu),kdx(w,mu),'b','LineWidth',2)

%plot(w,fragilityIDTFexample2(-0.9,-0.2,w,0.5),'b','LineWidth',2)
%axis([0 2 -.1 .2])

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

