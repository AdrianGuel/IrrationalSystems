kp=1;
kd=1;
t = 0:.01:30;
result1=euler_inversion(@(s) ((kp+kd*s^mu)/(sqrt(s^2-1)+(kp+kd*s^mu)))*(1/s), t);
plot(t,result1,'b','LineWidth',2)