a=2;omig=pi/4;R=3; t=0:0.1:10;
[yreal,ym]=funtrackingsnake(a,omig,t,R);
plot(t,yreal,t,ym)
xlabel('t')
ylabel('��������')
legend('Ŀ���˶���ʵ������','��������������')