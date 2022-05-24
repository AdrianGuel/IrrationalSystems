clear all; close all;

%mu=-0.8;
mu=0.4;
p=2;
q=2;
ml=1;
k=0.2E0;
be=0.4E0;
a=(p-1)*k;
b=(q-1)*be;
c=4*(p+q-1)*k*be;
m=2*ml;
w=0.0001:.01:2;
kp0=-5:1:3;
kd0=-ones(1,length(kp0));
%kp0(:)=-1;
%kpi=-3:1:4;
%kdi=-(m/(2*b))*ones(1,length(kpi));
 subplot(2,2,[1,3])
 kpm=-1:0.5:1;
 kdm=2*ones(1,length(kpm));
 H=area(kpm,kdm,-4);
 set(H,'FaceColor',[0.85 0.85 0.85],'LineStyle','none','ShowBaseLine','off');
 hold on
 w=0:.1:20;
 H3=area(kpx(w,mu,a,b,c,m),kdx(w,mu,a,b,c,m),-4);
 set(H3,'FaceColor',[1 1 1],'LineStyle','none','ShowBaseLine','off');
 hold on
%  w=0.0001:.1:2.5;
%  H3=area(kpx(w,mu,a,b,c,m),kdx(w,mu,a,b,c,m),-4);
%  set(H3,'FaceColor',[0.85 0.85 0.85],'LineStyle','none','ShowBaseLine','off');
%  hold on
%  H=area(kd0,kp0,-4);
%  set(H,'FaceColor',[1 1 1],'LineStyle','none','ShowBaseLine','off');
 w=0.0001:0.1:4;
plot(kpx(w,mu,a,b,c,m),kdx(w,mu,a,b,c,m),'b','LineWidth',2)
hold on
plot(kd0,kp0,'r','LineWidth',2)
% hold on
% plot(kpi,kdi,'y','LineWidth',2)
axis('square')
xlabel('$k_p(\omega)$','Interpreter','Latex','FontSize', 16)
ylabel('$k_\eta(\omega)$','Interpreter','Latex','FontSize', 16)
% axis([-1 0 -3 20]);
axis([-2 1 -3 2]);
 txt = 'Stability Region';
text(-0.3,15,txt,'FontSize',8)
 x = [.2 0.23];
 y = [.43 0.35];
 annotation('textarrow',x,y,'String','$RRB$','Interpreter','Latex','FontSize', 12)
 x = [.32 0.36];
 y = [.55 0.55];
 annotation('textarrow',x,y,'String','$CRB$','Interpreter','Latex','FontSize', 12)
 kp=[0.528 0 -1.5 0];
 kd=[-0.87 -2.3 -2 1];
 hold on
 plot(kp(1),kd(1),'kp','LineWidth',2)
 text(kp(1)-0.13,kd(1)+0.3,'$\mathbf{k}_1$','Interpreter','Latex','Color','k','FontSize',14)
hold on
 plot(kp(2),kd(2),'kp','LineWidth',2)
 text(kp(2)-0.13,kd(2)+0.3,'$\mathbf{k}_2$','Interpreter','Latex','Color','k','FontSize',14)
 hold on
 plot(kp(3),kd(3),'kp','LineWidth',2)
 text(kp(3)-0.13,kd(3)+0.3,'$\mathbf{k}_3$','Interpreter','Latex','Color','k','FontSize',14)
  hold on
 plot(kp(4),kd(4),'kp','LineWidth',2)
 text(kp(4)+0.05,kd(4)+0.1,'$\mathbf{k}_4$','Interpreter','Latex','Color','k','FontSize',14)
% 
w=0.0001:0.001:4;
subplot(2,2,2)
plot(kpx(w,mu,a,b,c,m),real(DSkp(w,mu,a,b,c,m,kpx(w,mu,a,b,c,m),kdx(w,mu,a,b,c,m))))
axis('square')
xlabel('$k_p(\omega)$','Interpreter','Latex','FontSize', 16)
ylabel('$\mathcal{S}_p$','Interpreter','Latex','FontSize', 16)
text(-0.8,-1,'$\omega$','Interpreter','Latex','Color','red','FontSize',14)
%axis([-1 -0.2 -1.5 -0.2]);
xlim([-2,1]);
subplot(2,2,4)
aux1=real(DSkd(w,mu,a,b,c,m,kpx(w,mu,a,b,c,m),kdx(w,mu,a,b,c,m)));
aux2=kdx(w,mu,a,b,c,m);
plot(kdx(w,mu,a,b,c,m),real(DSkd(w,mu,a,b,c,m,kpx(w,mu,a,b,c,m),kdx(w,mu,a,b,c,m))))
hold on
axis('square')
xlabel('$k_\eta(\omega)$','Interpreter','Latex','FontSize', 16)
ylabel('$\mathcal{S}_\eta$','Interpreter','Latex','FontSize', 16)
%[ d, ix ] = min(abs(aux1));
%plot(aux2(ix),aux1(ix),'x','LineWidth',3)
%axis([-.5 20 -0.5 0.5]);
xlim([-3,2])
text(1,-10,'$\omega$','Interpreter','Latex','Color','red','FontSize',14)
%text(aux2(ix),aux1(ix)-.1,'$\mathcal{S}_\eta=0$','Interpreter','Latex','Color','green','FontSize',10)
set(gcf,'color','w');

%pickup 4 points
%[kp,kd] =ginput(4);
% 
%Simulation impulse response
t = 0:.1:50;
result1=euler_inversion(@(s) ((kp(1)+kd(1)*s^mu)*(a+b*s+sqrt(c*s+(a+b*s)^2))/(m*s^2+b*s+a+sqrt((a+b*s)^2+c*s)+(kp(1)+kd(1)*s^mu)*(b*s+a+sqrt((a+b*s)^2+c*s))))*(1/s), t);
result2=euler_inversion(@(s) ((kp(2)+kd(2)*s^mu)*(a+b*s+sqrt(c*s+(a+b*s)^2))/(m*s^2+b*s+a+sqrt((a+b*s)^2+c*s)+(kp(2)+kd(2)*s^mu)*(b*s+a+sqrt((a+b*s)^2+c*s))))*(1/s), t);

figure
subplot(4,1,1)
plot(t,result1,'b','LineWidth',2)
%axis([-0.5 -0.2 -10 50]);
%axis('square')
%xlabel('Time (sec.)','Interpreter','Latex','FontSize', 16)
ylabel('$y(t)$','Interpreter','Latex','FontSize', 16)
str = sprintf('$k_1=(k_p,k_\\eta)=$ (%.2f,%.2f)',kp(1),kd(1));
title(str,'Interpreter','Latex')

subplot(4,1,2)
plot(t,result2,'b','LineWidth',2)
xlim([0,20])
%axis('square')
%xlabel('Time (sec.)','Interpreter','Latex','FontSize', 16)
ylabel('$y(t)$','Interpreter','Latex','FontSize', 16)
str = sprintf('$k_2=(k_p,k_\\eta)=$ (%.2f,%.2f)',kp(2),kd(2));
title(str,'Interpreter','Latex')
subplot(4,1,3)
t = 0:.1:50;
result3=euler_inversion(@(s) ((kp(3)+kd(3)*s^mu)*(a+b*s+sqrt(c*s+(a+b*s)^2))/(m*s^2+b*s+a+sqrt((a+b*s)^2+c*s)+(kp(3)+kd(3)*s^mu)*(b*s+a+sqrt((a+b*s)^2+c*s))))*(1/s), t);
plot(t,result3,'b','LineWidth',2)
xlim([0,11])
%axis('square')
%xlabel('Time (sec.)','Interpreter','Latex','FontSize', 16)
ylabel('$y(t)$','Interpreter','Latex','FontSize', 16)
str = sprintf('$k_3=(k_p,k_\\eta)=$ (%.2f,%.2f)',kp(3),kd(3));
title(str,'Interpreter','Latex')
subplot(4,1,4)
t = 0:.1:50;
result4=euler_inversion(@(s) ((kp(4)+kd(4)*s^mu)*(a+b*s+sqrt(c*s+(a+b*s)^2))/(m*s^2+b*s+a+sqrt((a+b*s)^2+c*s)+(kp(4)+kd(4)*s^mu)*(b*s+a+sqrt((a+b*s)^2+c*s))))*(1/s), t);
plot(t,result4,'b','LineWidth',2)
%axis('square')
xlabel('Time (sec.)','Interpreter','Latex','FontSize', 16)
ylabel('$y(t)$','Interpreter','Latex','FontSize', 16)
str = sprintf('$k_4=(k_p,k_\\eta)=$ (%.2f,%.2f)',kp(4),kd(4));
title(str,'Interpreter','Latex')
set(gcf,'color','w');

%fragility
%A=inftreefragility(1,1,2,a,b,c,m);

function kpaux=kpx(w,mu,a,b,c,m)
kpaux=(-1).*(b.^2.*(w.^2).^(1+(1/2).*mu).*sin(mu.*angle(sqrt(-1).*w))+a.^2.*( ...
  w.^2).^((1/2).*mu).*sin(mu.*angle(sqrt(-1).*w))+2.*a.*(w.^2).^((1/2).* ...
  mu).*((2.*a.*b.*w+c.*w).^2+(a.^2+(-1).*b.^2.*w.^2).^2).^(1/4).*cos((1/2) ...
  .*angle(sqrt(-1).*c.*w+(a+sqrt(-1).*b.*w).^2)).*sin(mu.*angle(sqrt(-1).* ...
  w))+(w.^2).^((1/2).*mu).*((2.*a.*b.*w+c.*w).^2+(a.^2+(-1).*b.^2.*w.^2) ...
  .^2).^(1/2).*cos((1/2).*angle(sqrt(-1).*c.*w+(a+sqrt(-1).*b.*w).^2)) ...
  .^2.*sin(mu.*angle(sqrt(-1).*w))+2.*b.*w.*(w.^2).^((1/2).*mu).*((2.*a.* ...
  b.*w+c.*w).^2+(a.^2+(-1).*b.^2.*w.^2).^2).^(1/4).*sin(mu.*angle(sqrt(-1) ...
  .*w)).*sin((1/2).*angle(sqrt(-1).*c.*w+(a+sqrt(-1).*b.*w).^2))+(w.^2).^( ...
  (1/2).*mu).*((2.*a.*b.*w+c.*w).^2+(a.^2+(-1).*b.^2.*w.^2).^2).^(1/2).* ...
  sin(mu.*angle(sqrt(-1).*w)).*sin((1/2).*angle(sqrt(-1).*c.*w+(a+sqrt(-1) ...
  .*b.*w).^2)).^2).^(-1).*((a+(-1).*m.*w.^2+((2.*a.*b.*w+c.*w).^2+(a.^2+( ...
  -1).*b.^2.*w.^2).^2).^(1/4).*cos((1/2).*angle(sqrt(-1).*c.*w+(a+sqrt(-1) ...
  .*b.*w).^2))).*(b.*w.*(w.^2).^((1/2).*mu).*cos(mu.*angle(sqrt(-1).*w))+ ...
  a.*(w.^2).^((1/2).*mu).*sin(mu.*angle(sqrt(-1).*w))+(w.^2).^((1/2).*mu) ...
  .*((2.*a.*b.*w+c.*w).^2+(a.^2+(-1).*b.^2.*w.^2).^2).^(1/4).*cos((1/2).* ...
  angle(sqrt(-1).*c.*w+(a+sqrt(-1).*b.*w).^2)).*sin(mu.*angle(sqrt(-1).*w) ...
  )+(w.^2).^((1/2).*mu).*((2.*a.*b.*w+c.*w).^2+(a.^2+(-1).*b.^2.*w.^2).^2) ...
  .^(1/4).*cos(mu.*angle(sqrt(-1).*w)).*sin((1/2).*angle(sqrt(-1).*c.*w+( ...
  a+sqrt(-1).*b.*w).^2)))+(-1).*(b.*w+((2.*a.*b.*w+c.*w).^2+(a.^2+(-1).* ...
  b.^2.*w.^2).^2).^(1/4).*sin((1/2).*angle(sqrt(-1).*c.*w+(a+sqrt(-1).*b.* ...
  w).^2))).*(a.*(w.^2).^((1/2).*mu).*cos(mu.*angle(sqrt(-1).*w))+(w.^2).^( ...
  (1/2).*mu).*((2.*a.*b.*w+c.*w).^2+(a.^2+(-1).*b.^2.*w.^2).^2).^(1/4).* ...
  cos(mu.*angle(sqrt(-1).*w)).*cos((1/2).*angle(sqrt(-1).*c.*w+(a+sqrt(-1) ...
  .*b.*w).^2))+(-1).*b.*w.*(w.^2).^((1/2).*mu).*sin(mu.*angle(sqrt(-1).*w) ...
  )+(-1).*(w.^2).^((1/2).*mu).*((2.*a.*b.*w+c.*w).^2+(a.^2+(-1).*b.^2.* ...
  w.^2).^2).^(1/4).*sin(mu.*angle(sqrt(-1).*w)).*sin((1/2).*angle(sqrt(-1) ...
  .*c.*w+(a+sqrt(-1).*b.*w).^2))));
end
function kdaux=kdx(w,mu,a,b,c,m)
kdaux=(-1).*m.*(w.^2).^(1+(-1/2).*mu).*csc(mu.*angle(sqrt(-1).*w)).*(b.*w+(( ...
  2.*a.*b.*w+c.*w).^2+(a.^2+(-1).*b.^2.*w.^2).^2).^(1/4).*sin((1/2).* ...
  angle(sqrt(-1).*c.*w+(a+sqrt(-1).*b.*w).^2))).*(a.^2+b.^2.*w.^2+2.*a.*(( ...
  2.*a.*b.*w+c.*w).^2+(a.^2+(-1).*b.^2.*w.^2).^2).^(1/4).*cos((1/2).* ...
  angle(sqrt(-1).*c.*w+(a+sqrt(-1).*b.*w).^2))+((2.*a.*b.*w+c.*w).^2+( ...
  a.^2+(-1).*b.^2.*w.^2).^2).^(1/2).*cos((1/2).*angle(sqrt(-1).*c.*w+(a+ ...
  sqrt(-1).*b.*w).^2)).^2+2.*b.*w.*((2.*a.*b.*w+c.*w).^2+(a.^2+(-1).* ...
  b.^2.*w.^2).^2).^(1/4).*sin((1/2).*angle(sqrt(-1).*c.*w+(a+sqrt(-1).*b.* ...
  w).^2))+((2.*a.*b.*w+c.*w).^2+(a.^2+(-1).*b.^2.*w.^2).^2).^(1/2).*sin(( ...
  1/2).*angle(sqrt(-1).*c.*w+(a+sqrt(-1).*b.*w).^2)).^2).^(-1);
end

function dskpaux=DSkp(w,mu,a,b,c,m,kp,kd)
dskpaux=((-1).*a+(sqrt(-1)*(-1)).*b.*w+(-1).*(sqrt(-1).*c.*w+(a+sqrt(-1).*b.*w) ...
  .^2).^(1/2)).*(b+(sqrt(-1)*2).*m.*w+(1/2).*(c+2.*b.*(a+sqrt(-1).*b.*w)) ...
  .*(sqrt(-1).*c.*w+(a+sqrt(-1).*b.*w).^2).^(-1/2)+(kp+kd.*(sqrt(-1).*w) ...
  .^mu).*(b+(1/2).*(c+2.*b.*(a+sqrt(-1).*b.*w)).*(sqrt(-1).*c.*w+(a+sqrt( ...
  -1).*b.*w).^2).^(-1/2))+kd.*mu.*(sqrt(-1).*w).^((-1)+mu).*(a+sqrt(-1).* ...
  b.*w+(sqrt(-1).*c.*w+(a+sqrt(-1).*b.*w).^2).^(1/2))).^(-1);
end

function dskdaux=DSkd(w,mu,a,b,c,m,kp,kd)
dskdaux=(-1).*(sqrt(-1).*w).^mu.*(a+sqrt(-1).*b.*w+(sqrt(-1).*c.*w+(a+sqrt(-1).* ...
  b.*w).^2).^(1/2)).*(b+(sqrt(-1)*2).*m.*w+(1/2).*(c+2.*b.*(a+sqrt(-1).* ...
  b.*w)).*(sqrt(-1).*c.*w+(a+sqrt(-1).*b.*w).^2).^(-1/2)+(kp+kd.*(sqrt(-1) ...
  .*w).^mu).*(b+(1/2).*(c+2.*b.*(a+sqrt(-1).*b.*w)).*(sqrt(-1).*c.*w+(a+ ...
  sqrt(-1).*b.*w).^2).^(-1/2))+kd.*mu.*(sqrt(-1).*w).^((-1)+mu).*(a+sqrt( ...
  -1).*b.*w+(sqrt(-1).*c.*w+(a+sqrt(-1).*b.*w).^2).^(1/2))).^(-1);
end

