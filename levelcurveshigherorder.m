clear all; close all; clc;


fighdl = figure(); view(3); hold on;
%view([0,0,1]);
set(fighdl,'defaulttextInterpreter','latex') 
set(gcf,'color','w');
axis square
xlabel('$k_p(\omega)$','FontSize', 16)
ylabel('$k_\eta(\omega)$','FontSize', 16)
zlabel('$\alpha$','FontSize', 16)
axis([-40 10 -5 30.05 -1 0]); 
kd0=0:1:10; w=0.00001:.01:8;
kp0=zeros(1,length(kd0));
kpi=-5:0.1:0.8;
kdi=kpi;
kdi(:)=-1;
z=ones(1,length(w)); z2=ones(1,length(kd0));

%mu=-0.5;
%plot(kpx(w,mu),kdx(w,mu),'--b',kd0,kp0,'r')
for mu=-1:.01:-0.1
[ d, ix ] = min( abs( kdx(w,mu)-30));
aux1=kdx(w,mu);
aux2=kpx(w,mu);
plot3(aux2,aux1,z.*mu,'b',kd0,kp0,z2.*mu,'r',[10 10],[0 30],[mu mu],'--b',[10 aux2(ix)],[30 aux1(ix)],[mu mu],'--b')
%plot(kpx(w,mu),kdx(w,mu),'b','LineWidth',3)
end
hold off

fighdl = figure(2); view(3); hold on;
%view([0,0,1]);
set(fighdl,'defaulttextInterpreter','latex') 
set(gcf,'color','w');
axis square
xlabel('$k_p(\omega)$','FontSize', 16)
ylabel('$k_\eta(\omega)$','FontSize', 16)
zlabel('$\alpha$','FontSize', 16)
axis([-1.3 2 -3 0.2 0 1]); 
w=0.00001:.01:80;
kpi=0.366:.1:2;
kdi=-ones(1,length(kpi));
z=ones(1,length(w)); 
z3=ones(1,length(kdi));
kd0=-2:.01:5; 
kp0=0.36602540378443865*ones(1,length(kd0));
z2=ones(1,length(kd0));
%mu=-0.5;
%plot(kpx(w,mu),kdx(w,mu),'--b',kd0,kp0,'r')
for mu=0.05:.05:.95
[ d, ix ] = min( abs( kpx(w,mu)-2));
aux1=kdx(w,mu);
aux2=kpx(w,mu);
if mu>=0.6
[x0,y0,segments]=selfintersect(aux2,aux1);
aux2r=aux2;
aux1r=aux1;
aux2r(segments(1):segments(2))=[];
aux1r(segments(1):segments(2))=[];
[ d2, ix2 ] = min( abs(kd0-aux1r(1)));
 if mu==0.8
plot3(aux2r,aux1r,ones(1,length(aux2r)).*mu,'b','LineWidth',3)
 plot3(kp0(ix2:end),kd0(ix2:end),z2(ix2:end).*mu,'r','LineWidth',3)
 plot3([0.3660 2],[0.2 0.2],[mu mu],'b','LineWidth',3)
 plot3([2 2],[0.2 aux1(ix)],[mu mu],'b','LineWidth',3)
 else
 plot3(aux2r,aux1r,ones(1,length(aux2r)).*mu,'b',kp0(ix2:end),kd0(ix2:end),z2(ix2:end).*mu,'r',[0.3660 2],[0.2 0.2],[mu mu],'--b',[2 2],[0.2 aux1(ix)],[mu mu],'--b')
 end
else
  [ d2, ix2 ] = min( abs(kd0-aux1(1)));
  plot3(aux2,aux1,ones(1,length(aux2)).*mu,'b',kp0(ix2:end),kd0(ix2:end),z2(ix2:end).*mu,'r',[0.3660 2],[0.2 0.2],[mu mu],'--b',[2 2],[0.2 aux1(ix)],[mu mu],'--b')
end
%plot(kpx(w,mu),kdx(w,mu),'b','LineWidth',3)
end
% X = [-1.5 -1 -1 -1.5 ];
% Y = [-1.5 -1.5  -1/sqrt(3) -1/sqrt(3)];
% Z = [0.5 0.5 0.5 0.5];
% fill3(X,Y,Z,[0.95 0.95 0.95])
kd0=-1:.1:0.2; 
kp0=0.36602540378443865*ones(1,length(kd0));
aux1=kdx(w,1);
aux2=kpx(w,1);
plot3(aux2,aux1,z.*1,'b','LineWidth',3)
plot3(kp0,kd0,ones(1,length(kd0)).*1,'r','LineWidth',3)
plot3(kpi,kdi,ones(1,length(kdi)).*1,'y','LineWidth',3)
plot3([0.3660 2],[0.2 0.2],[1 1],'b','LineWidth',3)
plot3([2 2],[0.2 -1],[1 1],'b','LineWidth',3)
plot3([-1.3 0.366],[-3 -3],[1 1],'k--','LineWidth',3)
plot3([-1.3 -1.3],[-3 -1],[1 1],'k--','LineWidth',3)
kd0=-3:.1:0.-1; 
kp0=0.36602540378443865*ones(1,length(kd0));
kpi=-1.3:.1:0;
kdi=-ones(1,length(kpi));
[ d2, ix2 ] = min( abs(kd0-aux1(1)));
[ d, ix ] = min(abs(kpi-aux2(end)));
plot3(kp0(1:ix2),kd0(1:ix2),ones(1,length(kd0(1:ix2))).*1,'r','LineWidth',3)
plot3(kpi(1:ix),kdi(1:ix),ones(1,length(kdi(1:ix))).*1,'y','LineWidth',3)
 hold off



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