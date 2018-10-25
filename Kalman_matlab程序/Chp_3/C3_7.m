clc
clear
z=[66.6,84.9,88.6,78.0,96.8,105.2,93.2,111.6,88.3,117.0,115.2]';
plot(z,'o--')
ylabel('�ֲ���������֣�')
k=1:11;
Hk=[k;ones(1,11)]';
estim=inv(Hk'*Hk)*Hk'*z;
ze=estim(1)*k+estim(2);

for i=1:4
ze(11+i)=estim(1)*(11+i)+estim(2)
 end

hold on
plot(ze,'y')

k=1:11;
Hk=[k.^2;k;ones(1,11)]';
estim=inv(Hk'*Hk)*Hk'*z;
ze=estim(1)*k.^2+estim(2)*k+estim(3);

 for i=1:4
ze(11+i)=estim(1)*((11+i).^2)+estim(2)*(11+i)+estim(3)
 end

hold on
plot(ze,'r')
k=1:11;
Hk=[k.^3;k.^2;k;ones(1,11)]';
estim=inv(Hk'*Hk)*Hk'*z;
ze=estim(1)*k.^3+estim(2)*k.^2+estim(3)*k+estim(4);

 for i=1:4
ze(11+i)=estim(1)*((11+i).^3)+estim(2)*((11+i).^2)+estim(3)*(11+i)+estim(4);
 end

hold on
plot(ze,'k')
k=1:11;
Hk=[k.^4;k.^3;k.^2;k;ones(1,11)]';
estim=inv(Hk'*Hk)*Hk'*z;
ze=estim(1)*k.^4+estim(2)*k.^3+estim(3)*k.^2+estim(4)*k+estim(5);
 
 for i=1:4
ze(11+i)=estim(1)*((11+i).^4)+estim(2)*((11+i).^3)+estim(3)*((11+i).^2)+estim(4)*(11+i)+estim(5);
 end

hold on
plot(ze,'g')
