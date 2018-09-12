clear
x_reality=[-0.7;1;1];           %��ʼ״̬
x_estimate=[0;0;0];             %��ʼ״̬�Ĺ���

Q=0.7;                          %����״̬Э����
R=1;                            % ��������Э����
P=[1 0 0;0 1 0;0 0 1];          %��ʼ���Ʒ���


n=3;            %ϵͳ��ά��
m=0.5;          %����ϵ��
L = 2 * n + 1;  %�ܵĲ�����ĸ���

for k=1:80
    x_reality(:,1)=[-0.7;1;1];
    x_reality(:,k+1)=[3*sin(2*x_reality(2,k));x_reality(1,k)+exp(-0.05*x_reality(3,k))+10;x_reality(1,k)*(x_reality(2,k)+x_reality(3,k))/5]+0.3*randn;
    x_array = x_reality;  %��ʵֵ����
    z(k)=x_reality(1,k)+x_reality(2,k)*x_reality(3,k)+0.5*randn;  %�۲�ֵ
end 

S=chol(P);
%%%%%%%%%״̬����
%%%%ѡ��ԳƲ���������״̬��sigma��
I=sqrt(n+m)*S;
x_sigma=x_estimate;

for i=2:n+1
    x_sigma(:,i) = x_estimate+I(:,i-1);
end

for i=n+2:L
    x_sigma(:,i) = x_estimate-I(:,i-n-1);
end

%%%%��Ӧ�ڸ�sigma���Ȩֵ Wm
w_1=m/(n+m);
w_2=1/(2*(n+m));


%%%%ʱ����¹���
%%sigama�㾭�������Ա任
for i=1:L
  r_sigma(:,i)=[3*sin(2*x_sigma(2,i));x_sigma(1,i)+exp(-0.05*x_sigma(3,i))+10;(x_sigma(1,i)*(x_sigma(2,i)+x_sigma(3,i)))/5+x_sigma(1,i)/2];
end

%%���ݺ�ľ�ֵ����ֵ��һ��Ԥ�⣩
x_next=zeros(3,1);
for i=1:L
    x_next=x_next+w_1*r_sigma(:,i);
end

%%���ݺ�ķ�������һ��Ԥ�⣩
 for i=1:L
    if i>1
        w=w_2;
    else 
        w=w_1;
    end
 end
 
 %%qr�ֽ� 
 for i=1:L
     A(:,i)=sqrt(w)*(r_sigma(:,i)-x_next);
 end
 
[Q S_next]=qr([A  sqrt(Q*eye(3))]',0);

%%cholupdate����
S_next=cholupdate(S_next,r_sigma(:,i)-x_next,'-') ;


%%%%������¹���
%%Sigma��ĸ���
I1=sqrt(n+m)*S_next;

 x_sigma_update=x_next;
 for i=2:n+1
     x_sigma_update(:,i)=x_next+I1(:,i-1);
 end
 for i=n+2:L
      x_sigma_update(:,i)=x_next-I1(:,i-n-1);
 end
 %%ͨ�����ⷽ��,�õ�һϵ�е�����
for i=1:n
     y_obser(i)=x_sigma_update(1,i)+x_sigma_update(2,i)*x_sigma_update(3,i);
end
 for i=n+1:L-1
     y_obser(i)=x_sigma_update(1,i-n)+x_sigma_update(2,i-n)*x_sigma_update(3,i-n);
 end
 for i=L;
     y_obser(i)=x_sigma_update(1,i-n-1)+x_sigma_update(2,i-n-1)*x_sigma_update(3,i-n-1);
 end
 %%�õ�һ��Ԥ�� 
z_next=0;
for i=1:L
     z_next=z_next+w_1*y_obser(i);
end
%%qr�ֽ�
for i=1:L
     B(:,i)=sqrt(w)*( y_obser(:,i)-z_next);
end
[R S_zz]=qr([B  sqrt(R)]',0);
%%cholupfdate����
S_zz=cholupdate(S_zz, y_obser(:,i)-z_next,'-');
 %%�õ�Э�������Pxz
 Pxz=zeros(3,1);
for i=1:L
   Pxz =Pxz+w*(x_sigma_update(:,i)- x_next)*(y_obser(i)- z_next)';
end
%%%%�����˲��������

%%����������
K = Pxz/(S_zz*S_zz');
U=K*S_zz;
%%�������ĸ���
for k=2:81
   x_estimate(:,k)=x_next+K*(z(k-1)-z_next);
end
x_estimate_array= [x_estimate];   %����ֵ����

%%����ĸ���
S=cholupdate(S_next,U,'-');

 %%������
 for k=1:80
     error(k)=0;
     for j=1:3
     error(k)=error(k)+(x_array(j,k) - x_estimate_array(j,k))^2;
     end
 end
 RMSE=sqrt(error/80);
 %%��ͼ

 figure
 k=1:1:80;
 plot(k,RMSE,'r-*');xlabel('���沽��');ylabel('RMSE');

   