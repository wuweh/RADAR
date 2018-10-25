function [xe,xee,pk,sss]=myUKFadfun(A,Q,R,xe,ym,mm,p,readerxy,sss)
%This function is to calculate the estimation state and the real state.
 
mx=6;
Xx=0.01;
Oo=2;
Vv=0;
   Kk=-3;
xe=[xe;zeros(size(xe));zeros(size(ym(1)))];
p=blkdiag(p,Q,R);
ssss=sqrtm((mx+Kk)*p); %%%先对方差矩阵求开方
if isreal(ssss) %%%%%%%%如果矩阵的方根为实数，则使用该数，
     sss=ssss;%%%否则，不更改由方差矩阵得到的数据
end
 
    s(:,1)=xe;
for ss=2:mx+1;
s(:,ss)=xe+sss(:,ss-1);
s(:,ss+mx)=xe-sss(:,ss-1);
end
 
Wx=ones(1,13)/2/(mx+Kk);
Wp=ones(1,13)/2/(mx+Kk);
 
Wx(1)=Kk/(mx+Kk);
Wp(1)=1-Xx*Xx+Oo+(Kk/(mx+Kk));
 
xee=0;
for i=1:length(s(1,:))
    xx(:,i)=A*s(1:6,i)+s(7:12,i); 
    xee=xee+Wx(i)*xx(:,i);
end
 
for i=1:length(s(1,:))
   xee1(:,i)=xx(:,i)-xee;
end
   
p=0;
for i=1:length(s(1,:))
   p=p+Wp(i)*xee1(:,i)*xee1(:,i)';
end
p=p+Q;
   
lmm=length(mm);   
 
for m=1:lmm
    for i=1:length(s(1,:))
        dmn(m,i)=sqrt((xx(1,i)-readerxy(1,mm(m))).^2+(xx(4,i)-readerxy(2,mm(m))).^2)+s(13,i);
    end
end
   
dm=0;
for i=1:length(s(1,:)) 
   dm=dm+Wx(i)*dmn(:,i);
end
   
for i=1:length(s(1,:))
   dm1(:,i)=dmn(:,i)-dm;
end
      
for m=1:lmm
   pxz=zeros(size(Wp(1)*xee1(:,1)*dm1(m,1)));
   for i=1:length(s(1,:))
        pxz=pxz+Wp(i)*xee1(:,i)*dm1(m,i)';
   end
   Pxz(:,m)=pxz;
end
   
for m1=1:lmm
    for m2=m1
       pzz=0;
       for i=1:length(s(1,:))
         pzz=pzz+Wp(i)*dm1(m1,i)*dm1(m2,i)';
       end
       Rm(m1,m2)=dm(m1)*dm(m2)*R;
       Pzz(m1,m2)=pzz+Rm(m1,m2);
     end
end            
K=Pxz*inv(Pzz);
pk=p-K*Pzz*K';
sumxee=0;
for  m=1:lmm      
    sumxee=sumxee+K(:,m)*(ym(m)-dm(m));
end
xe=xee+sumxee;  
