function [xe,xee,pk]=myEKFadfun(A,Q,R,xe,ym,mm,p,readerxy)
%This function is to calculate the estimation state and the real state.
   xee=A*xe;
   p=A*p*A'+Q;
  dm=sqrt((xee(1)-readerxy(1,mm)).^2+(xee(4)-readerxy(2,mm)).^2);
   Rm=dm.^2*R;
   for m=1:length(mm)
       Hm(m,:)=[(xee(1)-readerxy(1,mm(m)))/dm(m),0,0,(xee(4)-readerxy(2,mm(m)))/dm(m),0,0];
   end
   
   sumHm=0;
   for m=1:length(mm)
       sumHm=sumHm+Hm(m,:)'*inv(Rm(m))*Hm(m,:);
   end
   
   pk=inv(inv(p)+sumHm);
   
   Km=[];
   for m=1:length(mm)
       Km=[Km pk*Hm(m,:)'*inv(Rm(m))];
   end
   
   sumxe=0;
     for m=1:length(mm)
       sumxe=sumxe+Km(:,m)*(ym(m)-dm(m));
   end
   xe=xee+sumxe;

