function  [xx1,xxe1,P33,NN,qqxx,RR0x,RR0y]=funDataDrivenModelYWnonparell(TT,R,ax,qqx,ay,qqy,xe,p,y,N,readerxy)
xx1=[];qqxx=[];qq1=[];Ea1=[];RR0=[];RR1=[];P33=[];NN=[];
RR0x=[];RR0y=[];
xxe1=[];
Ea=[xe(3);xe(6)];
sss=ones(13);
MM=3;
[J,I]=size(y);
 
for i=1:I %￥￥￥￥
%%%%%%%%%%%%%%%%%%%%%%%%%%适用于不规则采样数据的周期
    if length(TT)==1
        T=TT;
    else
        T=TT(i);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 
ym=[];mm=[];
for j=1:J  %%找出第几个传感器有数据
   if isnan(y(j,i)) 
       
   else
       ym=[ym;y(j,i)];
       mm=[mm;j];
       end
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if i==1  %））））
        %%%%%%%%%%%%%从0到1步
[Ax,Qx,Ux]=myStarmodel(T,ax,qqx);
[Ay,Qy,Uy]=myStarmodel(T,ay,qqy);
A=blkdiag(Ax,Ay);
U=blkdiag(Ux,Uy);
Q=blkdiag(Qx,Qy);
 
 [xe,xee,p,sss]=myUKFadfun(A,Q,R,xe,ym,mm,p,readerxy,sss);%%%使用UKF方法进行非线性估计
 
%%%%%%%%%得到第1步的均值和两个自相关函数
Eax=xee(3);
R0x=xee(3)*xee(3);
R1x=xee(3);
 
Eay=xee(6);
R0y=xee(6)*xee(6);
R1y=xee(6);
 
Ea=[Eax;Eay];
    else%%%%其他步 %））））
%%%%%%%%%%%%%%%%%%%%%%%%% x轴
 
    if i>4   %* %%%%%第4步以后进行统计计算自相关函数，否则缺少必须的统计数据量
    bx=R1x/R0x; 
    if bx>0 
       qbx=(R0x-bx*R1x);
ax=-log(bx)/T;
 
qqx=2*ax*qbx/(1-bx*bx);
 
 if qqx>10000^2  qqx=10000^2;end  %%%%%%%%%%%%qq设置为有界值，如果qq发散，则系统发散
 end  %%如果b小于0，则会出现计算出复数的情况，出现计算问题，不更新模型参数a和qq
    end  %*    
%%%%%%%%%%%%%%%%%%%%%%%%%% y轴   
    if i>4   %* %%%%%第4步以后进行统计计算自相关函数，否则缺少必须的统计数据量
    by=R1y/R0y; 
    if by>0 
       qby=(R0y-by*R1y);
ay=-log(by)/T;
qqy=2*ay*qby/(1-by*by);
 if qqy>10000^2  qqy=10000^2;end  %%%%%%%%%%%%qq设置为有界值，如果qq发散，则系统发散
 end  %%如果b小于0，则会出现计算出复数的情况，出现计算问题，不更新模型参数a和qq
    %*   
    end  
    [Ax,Qx,Ux]=myStarmodel(T,ax,qqx);
[Ay,Qy,Uy]=myStarmodel(T,ay,qqy);
A=blkdiag(Ax,Ay);
U=blkdiag(Ux,Uy);
Q=blkdiag(Qx,Qy);
 
 
    lmm=length(mm);   
[xe,xee,p,sss]=myUKFadfun(A,Q,R,xe,ym,mm,p,readerxy,sss);
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%计算自适应模型  
     R0x=R0x+((xe(3))*(xe(3))-R0x)/i;    %%%%%%如果这样取，则利用估计值进行自适应参数计算
    R1x=R1x+((xe(3))*(xx1(3,i-1))-R1x)/i;  %%%%%不是并行的
    
    R0y=R0y+((xe(6))*(xe(6))-R0y)/i;
    R1y=R1y+((xe(6))*(xx1(6,i-1))-R1y)/i;

 end %））））
 

qqxx=[qqxx qqx];
RR0x=[RR0x R0x];
RR0y=[RR0y R0y];
 
xx1=[xx1 xe];
xxe1=[xxe1 xee];
P33=[P33 p(3,3)];
NN=[NN N];
end %￥￥￥￥
