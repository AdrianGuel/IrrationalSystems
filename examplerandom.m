t = 0:.001:20;
mu=0.5;
kp=10;
kd=10;
result1=euler_inversion(@(s) ((kp+kd*s^mu)/(sqrt(s-1)+(kp+kd*s^mu)))*(1/s), t);
plot(t,result1)