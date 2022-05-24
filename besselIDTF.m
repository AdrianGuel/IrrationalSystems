clear all;
close all;

mu=0.6;
w=0:.1:150;
kd0=-5:.1:10;
kp0=zeros(length(kd0));
kp0(:)=-1;

subplot(2,2,[1,3])
kpm=-1:1:1;
kdm=zeros(length(kpm));
kdm(:)=1;
H=area(kpm,kdm,-1);
set(H(1),'FaceColor',[0.85 0.85 0.85],'LineStyle','none','ShowBaseLine','off');
hold on
H=area(kpx(w,mu),kdx(w,mu),-1);
set(H(1),'FaceColor',[1 1 1],'LineStyle','none','ShowBaseLine','off');
hold on
plot(kpx(w,mu),kdx(w,mu),'b','LineWidth',3)
hold on
plot(kp0,kd0,'r','Linewidth',3)
txt = 'Stability region';
text(-0.7,0.85,txt)
x = [.33 0.37];
y = [.43 0.43];
annotation('textarrow',x,y,'String','$CRB$','Interpreter','Latex','FontSize', 12)
x = [.2 0.24];
y = [.7 0.7];
annotation('textarrow',x,y,'String','$RRB$','Interpreter','Latex','FontSize', 12)
x = [.4 0.45];
y = [.33 0.33];
annotation('textarrow',x,y,'Color','magenta')
text(0.2,-0.7,'$S_p<0$','Interpreter','Latex','Color','magenta','FontSize',11)
x = [.35 0.35];
y = [.5 0.57];
annotation('textarrow',x,y,'Color','magenta')
text(-0.8,.1,'$S_d<0$','Interpreter','Latex','Color','magenta','FontSize',11)
%plot(kpm,kdm);
kp=[0.3 -0.41 -1.69];
kd=[0.52 -0.73 0.34];
hold on
plot(kp(1),kd(1),'kp','LineWidth',2)
text(kp(1)-0.2,kd(1)+0.15,'$\mathbf{k}_1$','Interpreter','Latex','Color','k','FontSize',14)
hold on
plot(kp(2),kd(2),'kp','LineWidth',2)
text(kp(2)-0.2,kd(2)+0.15,'$\mathbf{k}_2$','Interpreter','Latex','Color','k','FontSize',14)
hold on
plot(kp(3),kd(3),'kp','LineWidth',2)
text(kp(3)-0.2,kd(3)+0.15,'$\mathbf{k}_3$','Interpreter','Latex','Color','k','FontSize',14)


axis('square')
xlabel('$k_p(\omega)$','Interpreter','Latex','FontSize', 16)
ylabel('$k_\eta(\omega)$','Interpreter','Latex','FontSize', 16)
 axis([-2 1 -1 1]);
% 
w=0.0001:0.1:15;
subplot(2,2,2)
plot(kpx(w,mu),real(DSkp(w,mu,kdx(w,mu))),'k--','LineWidth',1)
axis('square')
xlabel('$k_p(\omega)$','Interpreter','Latex','FontSize', 16)
ylabel('$\mathcal{S}_p$','Interpreter','Latex','FontSize', 16)
axis([-1 10 -1.4 0.2]);
x = [.75 0.8];
y = [.66 0.66];
annotation('textarrow',x,y,'Color','red')
text(5,-1,'$\omega$','Interpreter','Latex','Color','red','FontSize',14)
subplot(2,2,4)
plot(kdx(w,mu),real(DSkd(w,mu,kdx(w,mu))),'k--','LineWidth',1)
axis('square')
xlabel('$k_\eta(\omega)$','Interpreter','Latex','FontSize', 16)
ylabel('$\mathcal{S}_\eta$','Interpreter','Latex','FontSize', 16)
axis([-3.5 1 -30 1]);
x = [.8 0.8];
y = [.2 0.3];
annotation('textarrow',x,y,'Color','red')
text(0.2,-17,'$\omega$','Interpreter','Latex','Color','red','FontSize',14)
 set(gcf,'color','w');

% %pickup 4 points
% [kp,kd] = ginput(3);
%Simulation impulse response
t = 0:.01:25;
result1=euler_inversion(@(s) ((kp(1)+kd(1)*s^mu)/(sqrt(s^2+1)+(kp(1)+kd(1)*s^mu)))*(1/s), t);
result2=euler_inversion(@(s) ((kp(2)+kd(2)*s^mu)/(sqrt(s^2+1)+(kp(2)+kd(2)*s^mu)))*(1/s), t);
result3=euler_inversion(@(s) ((kp(3)+kd(3)*s^mu)/(sqrt(s^2+1)+(kp(3)+kd(3)*s^mu)))*(1/s), t);
%result4=euler_inversion(@(s) ((kp(4)+kd(4)*s^mu)/(sqrt(s+1)+(kp(4)+kd(4)*s^mu))), t);

figure
subplot(3,1,1)
plot(t,result1,'b','LineWidth',2)
%axis([-0.5 -0.2 -10 50]);
%axis('square')
%xlabel('Time (sec.)','Interpreter','Latex','FontSize', 14)
ylabel('$y(t)$','Interpreter','Latex','FontSize', 16)
str = sprintf('$k_1=(k_p,k_\\eta)$= (%.2f,%.2f)',kp(1),kd(1));
title(str,'Interpreter','Latex')

subplot(3,1,2)
plot(t,result2,'b','LineWidth',2)
%axis('square')
%xlabel('Time (sec.)','Interpreter','Latex','FontSize', 14)
ylabel('$y(t)$','Interpreter','Latex','FontSize', 16)
str = sprintf('$k_2=(k_p,k_\\eta)$= (%.2f,%.2f)',kp(2),kd(2));
title(str,'Interpreter','Latex')
subplot(3,1,3)
plot(t,result3,'b','LineWidth',2)
%axis('square')
xlabel('Time (sec.)','Interpreter','Latex','FontSize', 14)
ylabel('$y(t)$','Interpreter','Latex','FontSize', 16)
str = sprintf('$k_3=(k_p,k_\\eta)$= (%.2f,%.2f)',kp(3),kd(3));
title(str,'Interpreter','Latex')
%subplot(2,2,4)
%plot(t,result4)
set(gcf,'color','w');

function kpaux=kpx(w,mu)
kpaux=(-1).*(((-1)+w.^2).^2).^(1/4).*cos((1/2).*angle(1+(-1).*w.^2))+(((-1)+ ...
  w.^2).^2).^(1/4).*cot(mu.*angle(sqrt(-1).*w)).*sin((1/2).*angle(1+(-1).* ...
  w.^2));
end
function kdaux=kdx(w,mu)
kdaux=(-1).*(w.^2).^((-1/2).*mu).*(((-1)+w.^2).^2).^(1/4).*csc(mu.*angle(sqrt( ...
  -1).*w)).*sin((1/2).*angle(1+(-1).*w.^2));
end

function Dkpaux=DSkp(w,mu,kd)
Dkpaux=(-1).*(kd.*mu.*(sqrt(-1).*w).^((-1)+mu)+sqrt(-1).*w.*(1+(-1).*w.^2).^( ...
  -1/2)).^(-1);
end

function Dkdaux=DSkd(w,mu,kd)
Dkdaux=(-1).*(sqrt(-1).*w).^mu.*(kd.*mu.*(sqrt(-1).*w).^((-1)+mu)+sqrt(-1).*w.* ...
  (1+(-1).*w.^2).^(-1/2)).^(-1);
end


