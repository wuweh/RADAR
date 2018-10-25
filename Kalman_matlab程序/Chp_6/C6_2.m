a=2;omig=pi/4;R=3; t=0:0.1:10;
[yreal,ym]=funtrackingsnake(a,omig,t,R);
plot(t,yreal,t,ym)
xlabel('t')
ylabel('测量数据')
legend('目标运动的实际数据','传感器测量数据')