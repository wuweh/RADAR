clc
clear
close all
z=[66.6,84.9,88.6,78.0,96.8,105.2,93.2,111.6,88.3,117.0,115.2]';
plot(z,'o--')
ylabel('钢产量（百万吨）')
k=1:11;
Hk=[k.^2;k;ones(1,11)]';
estim=inv(Hk'*Hk)*Hk'*z
ze=estim(1)*k.^2+estim(2)*k+estim(3);

 for i=1:4
ze(11+i)=estim(1)*((11+i).^2)+estim(2)*(11+i)+estim(3);
 end

hold on
plot(ze)