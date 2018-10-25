%ʹ��˵����ʹ�����������Բ����켣�ĸ����㣬����Ҽ�Ϊ�����㡣
axis([0 100 0 100])
hold on
% ���������ʼ��
xy = [];
n = 0;
% ʹ��ѭ�����õ����������ǵ�����λ��
disp('Left mouse button picks points.')       %%������������ʾ��������
disp('Right mouse button picks last point.')
but = 1;
while but == 1
    [xi,yi,but] = ginput(1);         %%�ú������Եõ��������������
    plot(xi,yi,'ro')
    n = n+1;
    xy(:,n) = [xi;yi];
end
% ���ò�ֵ������ù⻬���ߣ�ģ��Ŀ���˶����ߺͲ������ݡ�
t = 1:n;
ts = 1: 0.1: n;
xys = spline(t,xy,ts);
plot(xys(1,:),xys(2,:),'b-'); %%����ɫֱ�߻���Ŀ���˶�����
Rx=10;Ry=10;    %%���ú��ᡢ����Ĳ�����������Rx��Ry��
plot(xys(1,:)+randn(size(xys(1,:)))*sqrt(Rx),xys(2,:)+randn(size(xys(1,:)))*sqrt(Ry),'k.');   %%����2ά�������ú�ɫ�㻭��ÿһ����������
xlabel('x��');ylabel('y��')
hold off
%����һ��ͼ�Ϸֱ𻭳��ᡢ�����Ŀ���˶����ݺͲ������ݡ�
figure
subplot(2,1,1),plot(xys(1,:))
hold on
subplot(2,1,1),plot(xys(1,:)+randn(size(xys(1,:)))*sqrt(Rx),'.')
ylabel('x��');
 
subplot(2,1,2),plot(xys(2,:))
hold on
subplot(2,1,2),plot(xys(2,:)+randn(size(xys(1,:)))*sqrt(Rx),'.')
ylabel('y��')
hold off
%�洢����
save  mytarget1 xys ts