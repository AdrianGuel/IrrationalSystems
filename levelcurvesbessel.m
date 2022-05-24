clear all; close all; clc;

fighdl = figure(); view(3); hold on;
%view([0,0,1]);
set(fighdl,'defaulttextInterpreter','latex') 
set(gcf,'color','w');
grid on
axis square
xlabel('$k_p(\omega)$','FontSize', 16)
ylabel('$k_\eta(\omega)$','FontSize', 16)
zlabel('$\alpha$','FontSize', 16)
axis([-5 10 -10 5 0 1.1]); 

kd0=0:.1:5; w=0:.01:100; 
kp0=-ones(1,length(kd0));
%kp0=zeros(1,length(kd0));
z=ones(1,length(w)); z2=ones(1,length(kd0));

%%

for mu=0.01:.1:1
[ d, ix ] = min( abs( kpx(w,mu)-10) );
aux1=kdx(w,mu);
plot3(kpx(w,mu),aux1,z.*mu,'--b',kp0,kd0,z2.*mu,'r',[-1 10],[5 5],[mu mu],'--b',[10 10],[5 aux1(ix)],[mu mu],'--b')
%plot(kpx(w,mu),kdx(w,mu),'b','LineWidth',3)
end
[ d, ix ] = min( abs( kpx(w,0.6)-10) );
aux1=kdx(w,0.6);
plot3(kpx(w,0.6),aux1,z.*0.6,'b','LineWidth',3)
plot3(kp0,kd0,z2.*0.6,'r','LineWidth',3)
plot3([-1 10],[5 5],[0.6 0.6],'b','LineWidth',3)
plot3([10 10],[5 aux1(ix)],[0.6 0.6],'b','LineWidth',3)
% X = [-5 -1 -1 -5 ];
% Y = [-10 -10  -1 -1];
% Z = [1 1 1 1];
% fill3(X,Y,Z,[0.95 0.95 0.95])
aux1=kdx(w,1);
aux2=kpx(w,1);
plot3(aux2,aux1,z.*1,'b','LineWidth',3)
plot3([-1 -1],[-10 -1],[1 1],'r','LineWidth',3)
plot3([-5 -1],[-1 -1],[1 1],'y','LineWidth',3)
plot3([-1 -1],[0 5],[1 1],'r','LineWidth',3)
plot3([aux2(end) 10],[-1 -1],[1 1],'y','LineWidth',3)
plot3([10 10],[0 5],[1 1],'--b')
plot3([10 -1],[5 5],[1 1],'--b')
plot3([-5 -1],[-10 -10],[1 1],'--k')
plot3([-5 -5],[-10 -1],[1 1],'--k')

%fill3([kpx(w,0.6),kp0,[-1 10],[10 10]],[aux1,kd0,[5 5],[5 aux1(ix)]],[z.*0.6,z2.*0.6,[0.6 0.6],[0.6 0.6]],[0.85 0.85 0.85])
hold off
%clearvars -except kd0  w

%%

% mu=0.001:.1:1;
% for j=1:length(mu),KPX(:,j)=kpx(w,mu(j)); KDX(:,j)=kdx(w,mu(j)); end
% ZM=z'*mu;
% 
% KP0=kp0'*ones(1,length(mu));
% KD0=kd0'*ones(1,length(mu));
% ZM2=z2'*mu;
% mesh(KPX,KDX,ZM, 'EdgeColor','b')
% mesh(KP0,KD0,ZM2,'EdgeColor','r')
%shading interp
%fill3([10*ones(size(KPX(:,end))); KPX(:,end)],[KDX(:,end); 5*ones(size(KDX(:,end)))],[ZM(:,end); ZM(:,end)],'b')
%hold off
% 
fighdl = figure(2); view(3); hold on;
%view([0,0,1]);
set(fighdl,'defaulttextInterpreter','latex') 
set(gcf,'color','w');
axis square
xlabel('$k_p(\omega)$','FontSize', 16)
ylabel('$k_\eta(\omega)$','FontSize', 16)
zlabel('$\alpha$','FontSize', 16)
axis([-10 5 -2 10 -1 0]); 

kd0=0:.1:5; w=1:.1:4000; 
kp0=zeros(1,length(kd0));
%kp0=zeros(1,length(kd0));
z=ones(1,length(w)); z2=ones(1,length(kd0));

%%

for mu=-1:.01:-.03
[ d, ix ] = min( abs( kdx(w,mu)-10) );
aux1=kdx(w,mu);
aux2=kpx(w,mu);
plot3(kpx(w,mu),aux1,z.*mu,'b',kd0,kp0,z2.*mu,'r',[5 5],[0 10],[mu mu],'--b',[5 aux2(ix)],[10 aux1(ix)],[mu mu],'--b')
%plot(kpx(w,mu),kdx(w,mu),'b','LineWidth',3)
end
axis([-11 5 -2 11 -1 0]); 
hold off


function kpaux=kpx(w,mu)
kpaux=(-1).*(((-1)+w.^2).^2).^(1/4).*cos((1/2).*angle(1+(-1).*w.^2))+(((-1)+ ...
  w.^2).^2).^(1/4).*cot(mu.*angle(sqrt(-1).*w)).*sin((1/2).*angle(1+(-1).* ...
  w.^2));
end
function kdaux=kdx(w,mu)
kdaux=(-1).*(w.^2).^((-1/2).*mu).*(((-1)+w.^2).^2).^(1/4).*csc(mu.*angle(sqrt( ...
  -1).*w)).*sin((1/2).*angle(1+(-1).*w.^2));
end