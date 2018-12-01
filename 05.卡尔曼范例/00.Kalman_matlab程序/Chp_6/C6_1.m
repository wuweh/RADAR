clc
clear
a=3;t=0:0.1:10;R=30;
[yreal,ym]=funtrackingline(a,t,R);
plot(t,yreal,t,ym)
xlabel('t')
ylabel('测量数据')
legend('目标运动的实际数据,传感器测量数据')