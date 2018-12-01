%对任意曲线进行跟踪
R=20;

%使用说明：使用鼠标左键可以产生轨迹的各个点，鼠标右键为结束点。
axis([0 100 0 100])
hold on
% Initially, the list of points is empty.
xy = [];
n = 0;
% Loop, picking up the points.
disp('Left mouse button picks points.')
disp('Right mouse button picks last point.')
but = 1;
while but == 1
    [xi,yi,but] = ginput(1);
    plot(xi,yi,'ro')
    n = n+1;
    xy(:,n) = [xi;yi];
end
% Interpolate with a spline curve and finer spacing.
t = 1:n;
ts = 1: 0.1: n;
xys = spline(t,xy,ts);
              
y=xys+sqrt(R)*randn(size(xys));
plot(xys(1,:),xys(2,:),'b-');
hold off
figure

%%估计横轴
xe=zeros(length(Q),1);p=10*eye(size(A));xx1=[];
for i=1:length(y(1,:))
[xe,p]=kalmanfun(A,C,Q,R,xe,y(1,i),p);
xx1=[xx1 xe];
end
%%%估计纵轴
xe=zeros(length(Q),1);p=10*eye(size(A));xx2=[];
for i=1:length(y(2,:))
[xe,p]=kalmanfun(A,C,Q,R,xe,y(2,i),p)
xx2=[xx2 xe];
end
plot(xys(1,:),xys(2,:),'b-');hold on
plot(y(1,:),y(2,:),'*');hold on,
plot(C*xx1,C*xx2,'o');hold off

diag(cov(xys'-[C*xx1;C*xx2]'))