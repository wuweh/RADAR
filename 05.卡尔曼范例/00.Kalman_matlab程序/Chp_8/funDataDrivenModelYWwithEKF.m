function  [xx1,xxe1,P33,NN]=funDataDrivenModelYWwithEKF(TT,R,ax,qqx,ay,qqy,xe,p,y,N,readerxy)
xx1=[];aa1=[];qq1=[];Ea1=[];RR0=[];RR1=[];P33=[];NN=[];
xxe1=[];
Ea=[xe(3);xe(6)];
 
MM=3;
[J,I]=size(y);
 
for i=1:I %�������� ѭ������ÿһ����������
%%%%%%%%%%%%%%%%%%%%%%%%%%���㲻����������ݵ�����
    if length(TT)==1
        T=TT;
    else
        T=TT(i);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�ҳ��ڼ����������в�������    
ym=[];mm=[];
for j=1:J  %%����ڼ����������в�������
   if isnan(y(j,i)) 
       
   else
       ym=[ym;y(j,i)];
       mm=[mm;j];
       end
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if i<=2  %%%%%%%%%%%%%%��0��1��
[Ax,Qx,Ux]=myStarmodel(T,ax,qqx);
[Ay,Qy,Uy]=myStarmodel(T,ay,qqy);
A=blkdiag(Ax,Ay);
U=blkdiag(Ux,Uy);
Q=blkdiag(Qx,Qy);
 
[xe,xee,p]=myEKFadfun(A,Q,R,xe,ym,mm,p,readerxy);
 
%%%%%%%%%�õ���1���ľ�ֵ����������غ���
Eax=xe(3);
R0x=xe(3)*xe(3);
R1x=xe(3);
 
Eay=xe(6);
R0y=xe(6)*xe(6);
R1y=xe(6);
 
Ea=[Eax;Eay];
    else%%%%������ %��������
%%%%%%%%%%%%%%%%%%%%%%%%% x��ļ���
  if i>10   %%%%%%��10���Ժ����ͳ�Ƽ�������غ���������ȱ�ٱ����ͳ��������
    bx=R1x/R0x; 
    if bx>0 
qbx=(R0x-bx*R1x);
ax=-log(bx)/T;
qqx=2*ax*qbx/(1-bx*bx);
 if qqx>10000^2  qqx=10000^2;end  %%%%%%%%%%%%qq����Ϊ�н�ֵ�����qq��ɢ����ϵͳ��ɢ
     end  %%���bС��0�������ּ������������������ּ������⣬������ģ�Ͳ���a��qq
  end  %*    
%%%%%%%%%%%%%%%%%%%%%%%%%% y��ļ���   
    if i>10  %* %%%%%��10���Ժ����ͳ�Ƽ�������غ���������ȱ�ٱ����ͳ��������
by=R1y/R0y; 
    if by>0 
       qby=(R0y-by*R1y);
ay=-log(by)/T;
qqy=2*ay*qby/(1-by*by);
    if qqy>10000^2  qqy=10000^2;end  %%%%%%%%%%%%qq����Ϊ�н�ֵ�����qq��ɢ����ϵͳ��ɢ
    end  %%���bС��0�������ּ������������������ּ������⣬������ģ�Ͳ���a��qq
     end %*   
    
    [Ax,Qx,Ux]=myStarmodel(T,ax,qqx);
[Ay,Qy,Uy]=myStarmodel(T,ay,qqy);
A=blkdiag(Ax,Ay);
U=blkdiag(Ux,Uy);
Q=blkdiag(Qx,Qy);
 
    [xe,xee,p]=myEKFadfun(A,Q,R,xe,ym,mm,p,readerxy); 

    R0x=R0x+((xe(3))*(xe(3))-R0x)/i;
    R1x=R1x+((xe(3))*(xx1(3,i-1))-R1x)/i;
    
    R0y=R0y+((xe(6))*(xe(6))-R0y)/i;
    R1y=R1y+((xe(6))*(xx1(6,i-1))-R1y)/i;
end %��������

xx1=[xx1 xe];
xxe1=[xxe1 xee];
P33=[P33 p(3,3)];
NN=[NN N];
end %��������
