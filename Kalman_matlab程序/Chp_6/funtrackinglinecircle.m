function [xt,y]=funtrackinglinecircle(T,Tt,Tz,R)
%T:采样周期
%Tt:直线运行时间
%Tz：圆周运行时间
%R：测量噪声方差
%本程序模拟匀速直线运动+圆周运动目标的真实轨迹
%运行描述：
% 针对某机动目标运动的特点，我们假设某运动目标经历了两个航行阶段：
% 初始匀速直线运动、匀速圆周转弯机动运动。将目标建立在二维坐标坐标系中，初始位置为XXX，目标线速度为XXX;
% 以XX度角向斜上方运动，目标运行41秒后，作向心wv ms/s^2的匀加速圆周运动，速度的大小v保持不变。
x0=[0 -2000]';%初始位置
v=200;%速度
vj=pi/4;%与水平方向所成角度

x1=x0;x=[];xt=[];
for t=1:T:Tt
    x2=[sin(vj)*v+x1(1);cos(vj)*v+x1(2)];
    x=[x x2];
    xt=[xt t];
    x1=x2;
end


% % figure% 圆周运动
x0=x2;
v=200;%线速度
T=1;
r=2000;%半径
w=v/r;
op=[x0(1)-r x0(2)];%圆心
x1=x0;www=[];ww=0;
for tt=t:Tz+t
ww=ww+w;
    x2=[cos(ww)*r+op(1);sin(ww)*r+op(2)];
    x=[x x2];
    xt=[xt tt];
    x1=x2;
%     www=[www ww];
end
y=x+randn(size(x))*sqrt(R);
