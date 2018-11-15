function [x]=IMM(A1,Q1,A2,Q2,A3,Q3,C,R,Pij,model0,y,x0)

%%%%%%%%%%%%%%%%%%%%%设置跟踪初值
 xei1=x0;pi1=10*eye(3);xxi1=[];
 xei2=x0;pi2=10*eye(3);xxi2=[];
 xei3=x0;pi3=10*eye(3);xxi3=[];
 x=[];pa=[];
 %%%%%%%%%%%%%%%%%%%%利用IMM进行跟踪
 [ii,jj]=size(Pij);

for i=1:length(y)
    for PPi=1:ii
 Pmodel(PPi,:)=Pij(PPi,:).*model0(PPi);
 end
 sumPm=sum(Pmodel);
 for PPj=1:jj
     model(:,PPj)=Pmodel(:,PPj)/sumPm(PPj);
 end
 
 
xej1=model(1,1)*xei1+model(2,1)*xei2+model(3,1)*xei3;
xej2=model(1,2)*xei1+model(2,2)*xei2+model(3,2)*xei3;
xej3=model(1,3)*xei1+model(2,3)*xei2+model(3,1)*xei3;

pj1=model(1,1)*(pi1+(xej1-xei1)*(xej1-xei1)')+model(2,1)*(pi2+(xej1-xei2)*(xej1-xei2)')+model(3,1)*(pi3+(xej1-xei3)*(xej1-xei3)');
pj2=model(1,2)*(pi1+(xej2-xei1)*(xej2-xei1)')+model(2,2)*(pi2+(xej2-xei2)*(xej2-xei2)')+model(3,2)*(pi3+(xej2-xei3)*(xej2-xei3)');
pj3=model(1,3)*(pi1+(xej3-xei1)*(xej3-xei1)')+model(2,3)*(pi2+(xej3-xei2)*(xej3-xei2)')+model(3,3)*(pi3+(xej3-xei3)*(xej3-xei3)');


[xej1,pkj1,v1,S1]=kalmanfunforIMM(A1,C,Q1,R,xej1,y(i),pj1);
[xej2,pkj2,v2,S2]=kalmanfunforIMM(A2,C,Q2,R,xej2,y(i),pj2);
[xej3,pkj3,v3,S3]=kalmanfunforIMM(A3,C,Q3,R,xej3,y(i),pj3);
% 
V1=exp(-v1*inv(S1)*v1)/sqrt(abs(2*pi*S1));
V2=exp(-v2*inv(S2)*v2)/sqrt(abs(2*pi*S2));
V3=exp(-v3*inv(S3)*v3)/sqrt(abs(2*pi*S3));

Vsum=sumPm*[V1,V2,V3]';

uj=sumPm.*[V1,V2,V3]/Vsum;

xe=xej1*uj(1)+xej2*uj(2)+xej3*uj(3);
p=pkj1+(xe-xej1)*(xe-xej1)'*uj(1)+pkj2+(xe-xej2)*(xe-xej2)'*uj(2)+pkj3+(xe-xej3)*(xe-xej3)'*uj(3);

model0=uj;
xei1=xej1;xei2=xej2;xei3=xej3;
pi1=pkj1;pi2=pkj2;pi3=pkj3;

x=[x xe];
xxi1=[xxi1 xei1];
xxi2=[xxi2 xei2];
xxi3=[xxi3 xei3];
 pa=[pa pkj1(1,1)];
end
% subplot(3,1,1),plot(y)
% subplot(3,1,2),plot(C*x)
% figure
% subplot(3,1,1),plot(C*xxi1)
% subplot(3,1,2),plot(C*xxi2)
% subplot(3,1,3),plot(C*xxi3)
% figure