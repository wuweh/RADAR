function [RR0,RR1,xx1,xxe1,aa1,qq1,Ea1]=StartrackingModel(y,a,qq,T,R,C)
 
xe=zeros(3,1);p=1000*eye(3);
xx1=[];aa1=[];qq1=[];Ea1=[];RR0=[];RR1=[];
xxe1=[];
for i=1:1  %%%%%%%%%%%%%从第0第1步
Ea=xe(3);
[A,Q,U]=myStarmodel(T,a,qq);
[xe,xee,p]=mykalmanadfun(A,U,C,Q,R,xe,y(i),p,Ea);
xx1=[xx1 xe];
xxe1=[xxe1 xee];
qsum=0;
%得到第1步的均值和两个自相关函数
Ea=xee(3);
R0=xee(3)*xee(3);
R1=xee(3)*0;
end
 
for i=2:length(y)
    if i>4
if R1/R0>0
    b=R1/R0;
      
     qb=(R0-b*R1);
a=-log(b)/T;
if 2*a*qb/(1-b*b)<100000;
qq=2*a*qb/(1-b*b);
    end
end
    end
    
    [A,Q,U]=myStarmodel(T,a,qq);
[xe,xee,p]=mykalmanadfun(A,U,C,Q,R,xe,y(i),p,Ea);
  Ea=((i-1)/i)*Ea+xee(3)/i;
 R0=R0+(10*xee(3)*xee(3)-R0)/i;
    R1=R1+(10*xee(3)*xxe1(3,i-1)-R1)/i;

RR0=[RR0 R0];
RR1=[RR1 R1];
xx1=[xx1 xe];
xxe1=[xxe1 xee];
aa1=[aa1 a];
qq1=[qq1 qq];
Ea1=[Ea1 Ea];
end