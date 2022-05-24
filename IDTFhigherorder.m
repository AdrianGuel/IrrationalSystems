mu=0.8;

kd0=-1.8:.1:-1;
kp0=kd0;
kp0(:)=0.36602540378443865;

subplot(2,2,[1,3])
kpm=0.366:0.1:1.7;
kdm=zeros(length(kpm));
kdm(:)=-1;
H=area(kpm,kdm,-3);
set(H,'FaceColor',[0.85 0.85 0.85],'LineStyle','none','ShowBaseLine','off');
hold on
w=0.00001:.1:0.48;
H3=area(kpx(w,mu),kdx(w,mu),-3);
set(H3,'FaceColor',[1 1 1],'LineStyle','none','ShowBaseLine','off');
hold on
w=3.2:.1:10;
H3=area(kpx(w,mu),kdx(w,mu),-3);
set(H3,'FaceColor',[1 1 1],'LineStyle','none','ShowBaseLine','off');
hold on
w=0.00001:.001:10;
plot(kpx(w,mu),kdx(w,mu),'b','LineWidth',2)
hold on
plot(kp0,kd0,'r','LineWidth',2)
axis('square')
xlabel('$k_p(\omega)$','Interpreter','Latex','FontSize', 16)
ylabel('$k_\eta(\omega)$','Interpreter','Latex','FontSize', 16)
axis([-0.3 1.5 -1.8 -1]);
%axis([-0.8 0.5 -2 5]);

x = [.35 0.28];
y = [.5 0.5];
annotation('textarrow',x,y,'String','$CRB$','Interpreter','Latex','FontSize', 12)
x = [.2 0.24];
y = [.7 0.7];
annotation('textarrow',x,y,'String','$RRB$','Interpreter','Latex','FontSize', 12)

%pickup 4 points
 %[kp,kd] = ginput(6);
kp=[0.95 0.94  -0.01 0.27 0.44 0.40];
kd=[-1.27 -1.66 -1.36 -1.55 -1.69 -1.38];
hold on
plot(kp(1),kd(1),'kp','LineWidth',2)
text(kp(1)+0.03,kd(1)+0.03,'$\mathbf{k}_1$','Interpreter','Latex','Color','k','FontSize',14)
hold on
plot(kp(2),kd(2),'kp','LineWidth',2)
text(kp(2)+0.03,kd(2)+0.03,'$\mathbf{k}_2$','Interpreter','Latex','Color','k','FontSize',14)
hold on
plot(kp(3),kd(3),'kp','LineWidth',2)
text(kp(3)+0.03,kd(3)+0.03,'$\mathbf{k}_3$','Interpreter','Latex','Color','k','FontSize',14)
hold on
plot(kp(4),kd(4),'kp','LineWidth',2)
text(kp(4)+0.03,kd(4)+0.03,'$\mathbf{k}_4$','Interpreter','Latex','Color','k','FontSize',14)
hold on
plot(kp(5),kd(5),'kp','LineWidth',2)
text(kp(5)+0.03,kd(5)+0.03,'$\mathbf{k}_5$','Interpreter','Latex','Color','k','FontSize',14)
hold on
plot(kp(6),kd(6),'kp','LineWidth',2)
text(kp(6)+0.03,kd(6)+0.03,'$\mathbf{k}_6$','Interpreter','Latex','Color','k','FontSize',14)
txt = 'Stability Region';
text(-1,-5,txt,'FontSize',12)

subplot(2,2,2)
plot(kpx(w,mu),real(DSkp(w,mu,kpx(w,mu),kdx(w,mu))),'k--')
axis('square')
xlabel('$k_p(\omega)$','Interpreter','Latex','FontSize', 16)
ylabel('$\mathcal{S}_p$','Interpreter','Latex','FontSize', 16)
axis([0 1.5 -4 3]);
x = [.75 0.7];
y = [.4 0.4];
annotation('textarrow',x,y,'Color','red')


subplot(2,2,4)
plot(kdx(w,mu),real(DSkd(w,mu,kpx(w,mu),kdx(w,mu))),'k--')
axis('square')
xlabel('$k_\eta(\omega)$','Interpreter','Latex','FontSize', 16)
ylabel('$\mathcal{S}_\eta$','Interpreter','Latex','FontSize', 16)
axis([-2 -0.3 -33 3]);
x = [.67 0.7];
y = [.75 0.65];
annotation('textarrow',x,y,'Color','red')
x = [.7 0.68];
y = [.76 0.85];
annotation('textarrow',x,y,'Color','red')
x = [.67 0.65];
y = [.3 0.2];
annotation('textarrow',x,y,'Color','red')
text(-0.9,-2,'$\omega$','Interpreter','Latex','Color','red','FontSize',14)
text(-1.4,45,'$\omega$','Interpreter','Latex','Color','red','FontSize',14)
set(gcf,'color','w');



% 
t = 0:.1:50;
result1=euler_inversion(@(s) ((s^2+2*s+1+sqrt(2*s+3))*(kp(1)+kd(1)*s^mu)/((s^3+3*s^2+4*s-2+sqrt(s+1))+(kp(1)+kd(1)*s^mu)*(s^2+2*s+1+sqrt(2*s+3))))*(1/s), t);
result2=euler_inversion(@(s) ((s^2+2*s+1+sqrt(2*s+3))*(kp(2)+kd(2)*s^mu)/((s^3+3*s^2+4*s-2+sqrt(s+1))+(kp(2)+kd(2)*s^mu)*(s^2+2*s+1+sqrt(2*s+3))))*(1/s), t);
result3=euler_inversion(@(s) ((s^2+2*s+1+sqrt(2*s+3))*(kp(3)+kd(3)*s^mu)/((s^3+3*s^2+4*s-2+sqrt(s+1))+(kp(3)+kd(3)*s^mu)*(s^2+2*s+1+sqrt(2*s+3))))*(1/s), t);
result4=euler_inversion(@(s) ((s^2+2*s+1+sqrt(2*s+3))*(kp(4)+kd(4)*s^mu)/((s^3+3*s^2+4*s-2+sqrt(s+1))+(kp(4)+kd(4)*s^mu)*(s^2+2*s+1+sqrt(2*s+3))))*(1/s), t);
result5=euler_inversion(@(s) ((s^2+2*s+1+sqrt(2*s+3))*(kp(5)+kd(5)*s^mu)/((s^3+3*s^2+4*s-2+sqrt(s+1))+(kp(5)+kd(5)*s^mu)*(s^2+2*s+1+sqrt(2*s+3))))*(1/s), t);
result6=euler_inversion(@(s) ((s^2+2*s+1+sqrt(2*s+3))*(kp(6)+kd(6)*s^mu)/((s^3+3*s^2+4*s-2+sqrt(s+1))+(kp(6)+kd(6)*s^mu)*(s^2+2*s+1+sqrt(2*s+3))))*(1/s), t);
figure
subplot(3,2,1)
plot(t,result1,'b','LineWidth',2)
%axis([-0.5 -0.2 -10 50]);
%axis('square')
%xlabel('Time (sec.)','Interpreter','Latex','FontSize', 16)
ylabel('$y(t)$','Interpreter','Latex','FontSize', 12)
str = sprintf('${k}_1=(k_p,k_\\eta)=$ (%.2f,%.2f)',kp(1),kd(1));
title(str,'Interpreter','Latex','position',[17 2]);
axis([0 40 -3 2]);

subplot(3,2,2)
plot(t,result2,'b','LineWidth',2)
%axis('square')
%xlabel('Time (sec.)','Interpreter','Latex','FontSize', 16)
ylabel('$y(t)$','Interpreter','Latex','FontSize', 12)
str = sprintf('${k}_2=(k_p,k_\\eta)=$ (%.2f,%.2f)',kp(2),kd(2));
title(str,'Interpreter','Latex','position',[5 2]);
axis([0 10.2 -3e11 2e11]);

subplot(3,2,3)
plot(t,result3,'b','LineWidth',2)
%axis('square')
%xlabel('Time (sec.)','Interpreter','Latex','FontSize', 16)
ylabel('$y(t)$','Interpreter','Latex','FontSize', 12)
str = sprintf('${k}_3=(k_p,k_\\eta)=$ (%.2f,%.2f)',kp(3),kd(3));
title(str,'Interpreter','Latex','position',[8 1]);
axis([0 12 -1e11 1e5]);

subplot(3,2,4)
plot(t,result4,'b','LineWidth',2)
%axis('square')
%xlabel('Time (sec.)','Interpreter','Latex','FontSize', 16)
ylabel('$y(t)$','Interpreter','Latex','FontSize', 12)
str = sprintf('${k}_4=(k_p,k_\\eta)=$ (%.2f,%.2f)',kp(4),kd(4));
title(str,'Interpreter','Latex','position',[5 1]);
axis([0 8 -1e11 1e5]);

subplot(3,2,5)
plot(t,result5,'b','LineWidth',2)
%axis('square')
xlabel('Time (sec.)','Interpreter','Latex','FontSize', 16)
ylabel('$y(t)$','Interpreter','Latex','FontSize', 12)
str = sprintf('${k}_5=(k_p,k_\\eta)=$ (%.2f,%.2f)',kp(5),kd(5));
title(str,'Interpreter','Latex','position',[2.2 1]);
axis([0 3.5 -1e11 1e5]);

subplot(3,2,6)
plot(t,result6,'b','LineWidth',2)
%axis('square')
xlabel('Time (sec.)','Interpreter','Latex','FontSize', 14)
ylabel('$y(t)$','Interpreter','Latex','FontSize', 12)
str = sprintf('${k}_5=(k_p,k_\\eta)=$ (%.2f,%.2f)',kp(6),kd(6));
title(str,'Interpreter','Latex','position',[25 545]);
%axis([0 40 -1 500]);
set(gcf,'color','w');

function kpaux=kpx(w,mu)
kpaux=(-1).*csc((1/2).*mu.*pi).*(1+2.*w.^2+w.^4+(9+4.*w.^2).^(1/2)+2.*(9+4.* ...
  w.^2).^(1/4).*((1+(-1).*w.^2).*cos((1/2).*atan((2/3).*w))+2.*w.*sin(( ...
  1/2).*atan((2/3).*w)))).^(-1).*(((-2)+(-3).*w.^2+(1+w.^2).^(1/4).*cos(( ...
  1/2).*atan(w))).*(2.*w.*cos((1/2).*mu.*pi)+(-1).*((-1)+w.^2).*sin((1/2) ...
  .*mu.*pi)+(9+4.*w.^2).^(1/4).*sin((1/2).*(mu.*pi+atan((2/3).*w))))+(-1) ...
  .*((1+(-1).*w.^2).*cos((1/2).*mu.*pi)+(9+4.*w.^2).^(1/4).*cos((1/2).*( ...
  mu.*pi+atan((2/3).*w)))+(-2).*w.*sin((1/2).*mu.*pi)).*(4.*w+(-1).*w.^3+( ...
  1+w.^2).^(1/4).*sin((1/2).*atan(w))));
end
function kdaux=kdx(w,mu)
kdaux=(-1).*w.^((-1).*mu).*csc((1/2).*mu.*pi).*(1+2.*w.^2+w.^4+(9+4.*w.^2).^( ...
  1/2)+2.*(9+4.*w.^2).^(1/4).*((1+(-1).*w.^2).*cos((1/2).*atan((2/3).*w))+ ...
  2.*w.*sin((1/2).*atan((2/3).*w)))).^(-1).*(8.*w+w.^3+w.^5+(-1).*w.*((-4) ...
  +w.^2).*(9+4.*w.^2).^(1/4).*cos((1/2).*atan((2/3).*w))+(-2).*w.*(1+w.^2) ...
  .^(1/4).*cos((1/2).*atan(w))+2.*(9+4.*w.^2).^(1/4).*sin((1/2).*atan(( ...
  2/3).*w))+3.*w.^2.*(9+4.*w.^2).^(1/4).*sin((1/2).*atan((2/3).*w))+(-1).* ...
  (1+w.^2).^(1/4).*(9+4.*w.^2).^(1/4).*sin((1/2).*(atan((2/3).*w)+(-1).* ...
  atan(w)))+(1+w.^2).^(1/4).*sin((1/2).*atan(w))+(-1).*w.^2.*(1+w.^2).^( ...
  1/4).*sin((1/2).*atan(w)));
end
function dskpaux=DSkp(w,mu,kp,kd)
dskpaux=(4+(1/2).*(1+sqrt(-1).*w).^(-1/2)+kp.*(2+(3+(sqrt(-1)*2).*w).^(-1/2)+( ...
  sqrt(-1)*2).*w)+kd.*(2+mu.*(2+(sqrt(-1)*(-1)).*(1+(3+(sqrt(-1)*2).*w).^( ...
  1/2)).*w.^(-1)+sqrt(-1).*w)+(3+(sqrt(-1)*2).*w).^(-1/2)+(sqrt(-1)*2).*w) ...
  .*(sqrt(-1).*w).^mu+(sqrt(-1)*6).*w+(-3).*w.^2).^(-1).*((-1)+(-1).*(3+( ...
  sqrt(-1)*2).*w).^(1/2)+w.*((sqrt(-1)*(-2))+w));
end
function dskdaux=DSkd(w,mu,kp,kd)
dskdaux=(-1).*(sqrt(-1).*w).^mu.*(4+(1/2).*(1+sqrt(-1).*w).^(-1/2)+kp.*(2+(3+( ...
  sqrt(-1)*2).*w).^(-1/2)+(sqrt(-1)*2).*w)+kd.*(2+mu.*(2+(sqrt(-1)*(-1)).* ...
  (1+(3+(sqrt(-1)*2).*w).^(1/2)).*w.^(-1)+sqrt(-1).*w)+(3+(sqrt(-1)*2).*w) ...
  .^(-1/2)+(sqrt(-1)*2).*w).*(sqrt(-1).*w).^mu+(sqrt(-1)*6).*w+(-3).*w.^2) ...
  .^(-1).*(1+(3+(sqrt(-1)*2).*w).^(1/2)+(-1).*w.*((sqrt(-1)*(-2))+w));
end