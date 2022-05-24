clear all; close all;
p=2;
q=2;
ml=1;
k=0.2E0;
be=0.4E0;
a=(p-1)*k;
b=(q-1)*be;
c=4*(p+q-1)*k*be;
m=2*ml;
kd0=-1:1:5;
kp0=zeros(1,length(kd0));
%kp0(:)=-1;
%kpi=-3:1:4;
%kdi=-(m/(2*b))*ones(1,length(kpi));


fighdl = figure(); view(3); hold on;
%view([0,0,1]);
set(fighdl,'defaulttextInterpreter','latex') 
set(gcf,'color','w');
axis square
xlabel('$k_p(\omega)$','FontSize', 16)
ylabel('$k_\eta(\omega)$','FontSize', 16)
zlabel('$\alpha$','FontSize', 16)
axis([-20 3 -3 20.05 -1 0]); 

w=0.0001:.01:5;
%kp0=zeros(1,length(kd0));
z=ones(1,length(w)); z2=ones(1,length(kd0));

%%

for mu=-1:.1:-.1
[ d, ix ] = min( abs(kdx(w,mu,a,b,c,m)-20) );
aux1=kdx(w,mu,a,b,c,m);
aux2=kpx(w,mu,a,b,c,m);
plot3(aux2,aux1,z.*mu,'b',kd0,kp0,z2.*mu,'r',[3 3],[0 20],[mu mu],'--b',[3 aux2(ix)],[20 aux1(ix)],[mu mu],'--b')
%plot(kpx(w,mu),kdx(w,mu),'b','LineWidth',3)
end
% mu=-0.8;
% [ d, ix ] = min( abs(kdx(w,mu,a,b,c,m)-20) );
% aux1=kdx(w,mu,a,b,c,m);
% aux2=kpx(w,mu,a,b,c,m);
%  plot3(aux2,aux1,z.*mu,'b','LineWidth',3)
%  plot3(kd0,kp0,z2.*mu,'r','LineWidth',3)
%  plot3([3 aux2(ix)],[20 aux1(ix)],[mu mu],'b','LineWidth',3)
%  plot3([3 3],[0 20],[mu mu],'b','LineWidth',3)
hold off
% 

fighdl = figure(2); view(3); hold on;
%view([0,0,1]);
set(fighdl,'defaulttextInterpreter','latex') 
set(gcf,'color','w');
%grid on
xlabel('$k_p(\omega)$','FontSize', 16)
ylabel('$k_\eta(\omega)$','FontSize', 16)
zlabel('$\alpha$','FontSize', 16)
%axis([-5 10 -15 5 0 1.1]); 

axis equal;
axis square;
axis([-1.4, 1.5,-3.5, 1.5, 0 ,1.1]);

w=0.0001:.01:20;
kd0=0:1:5;
kp0=-ones(1,length(kd0));
z=ones(1,length(w)); z2=ones(1,length(kd0));
kpi=-15:1:10;
kdi=-(m/(2*b))*ones(1,length(kpi));
z3=ones(1,length(kpi));
%%

for mu=0.01:.05:0.99
[ d, ix ] = min( abs(kpx(w,mu,a,b,c,m)-1.5) );
aux1=kdx(w,mu,a,b,c,m);
aux2=kpx(w,mu,a,b,c,m);
%plot3(aux2,aux1,z.*mu,'b',kp0,kd0,z2.*mu,'r',[-1 3],[1.5 1.5],[mu mu],'b-',[3 3],[1.5 aux1(ix)],[mu mu],'b-')
p1 = plot3(aux2,aux1,z.*mu,'b-','LineWidth',1); 
p1.Color(4) = 0.25;
p2 = plot3(kp0,kd0,z2.*mu,'r-','LineWidth',1); 
p2.Color(4) = 0.25;
p3 = plot3([-1 3],[1.5 1.5],[mu mu],'b-','LineWidth',1); 
p3.Color(4) = 0.25;
p4 = plot3([1.5 1.5],[1.5 aux1(ix)],[mu mu],'b-','LineWidth',1); 
p4.Color(4) = 0.25;
%plot(kpx(w,mu),kdx(w,mu),'b','LineWidth',3)
end
%X = [-5 -1 -1 -5 ];
%Y = [-15 -15  -(m/(2*b)) -(m/(2*b))];
%Z = [1 1 1 1];
%fill3(X,Y,Z,[0.95 0.95 0.95])
 mu=1;
 kd0=-15:1:5;
kp0=-ones(1,length(kd0));
z2=ones(1,length(kd0));
% [ d, ix ] = min( abs(kpx(w,mu,a,b,c,m)-20) );
 aux1=kdx(w,mu,a,b,c,m);
 aux2=kpx(w,mu,a,b,c,m);
  plot3(aux2,aux1,z.*mu,'b','LineWidth',3)
  plot3([-1 -1],[-15 -(m/(2*b))],[mu mu],'r','LineWidth',3)
   plot3([-1 -1],[0 5],[mu mu],'r','LineWidth',3)
 plot3([-5 -1],[-(m/(2*b)) -(m/(2*b))],[mu mu],'y','LineWidth',3)
plot3([aux2(end) 10],[-(m/(2*b)) -(m/(2*b))],[mu mu],'y','LineWidth',3)
  plot3([-1 10],[5 5],[mu mu],'b','LineWidth',3)
  plot3([10 10],[5 -(m/(2*b))],[mu mu],'b','LineWidth',3)
 plot3([-1.4 -1],[-3.5 -3.5],[mu mu],'--k','LineWidth',2)
 plot3([-1.4 -1.4],[-3.5 -(m/(2*b))],[mu mu],'--k')

hold on
r = 0.3493235542;
[X,Y,Z] = sphere;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;

h=surf(X2+0.528,Y2-0.87,Z2+0.4)
set(h, 'FaceAlpha', 0.5)
%shading interp

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
