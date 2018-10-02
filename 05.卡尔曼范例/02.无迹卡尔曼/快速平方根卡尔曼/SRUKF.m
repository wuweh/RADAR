clear
x_reality=[-0.7;1;1];           %初始状态
x_estimate=[0;0;0];             %初始状态的估计

Q=0.7;                          %过程状态协方差
R=1;                            % 测量噪声协方差
P=[1 0 0;0 1 0;0 0 1];          %初始估计方差


n=3;            %系统的维数
m=0.5;          %比例系数
L = 2 * n + 1;  %总的采样点的个数

for k=1:80
    x_reality(:,1)=[-0.7;1;1];
    x_reality(:,k+1)=[3*sin(2*x_reality(2,k));x_reality(1,k)+exp(-0.05*x_reality(3,k))+10;x_reality(1,k)*(x_reality(2,k)+x_reality(3,k))/5]+0.3*randn;
    x_array = x_reality;  %真实值数组
    z(k)=x_reality(1,k)+x_reality(2,k)*x_reality(3,k)+0.5*randn;  %观测值
end 

S=chol(P);
%%%%%%%%%状态估计
%%%%选择对称采样，构造状态的sigma点
I=sqrt(n+m)*S;
x_sigma=x_estimate;

for i=2:n+1
    x_sigma(:,i) = x_estimate+I(:,i-1);
end

for i=n+2:L
    x_sigma(:,i) = x_estimate-I(:,i-n-1);
end

%%%%对应于各sigma点的权值 Wm
w_1=m/(n+m);
w_2=1/(2*(n+m));


%%%%时间更新过程
%%sigama点经过非线性变换
for i=1:L
  r_sigma(:,i)=[3*sin(2*x_sigma(2,i));x_sigma(1,i)+exp(-0.05*x_sigma(3,i))+10;(x_sigma(1,i)*(x_sigma(2,i)+x_sigma(3,i)))/5+x_sigma(1,i)/2];
end

%%传递后的均值（均值的一步预测）
x_next=zeros(3,1);
for i=1:L
    x_next=x_next+w_1*r_sigma(:,i);
end

%%传递后的方差（方差的一步预测）
 for i=1:L
    if i>1
        w=w_2;
    else 
        w=w_1;
    end
 end
 
 %%qr分解 
 for i=1:L
     A(:,i)=sqrt(w)*(r_sigma(:,i)-x_next);
 end
 
[Q S_next]=qr([A  sqrt(Q*eye(3))]',0);

%%cholupdate更新
S_next=cholupdate(S_next,r_sigma(:,i)-x_next,'-') ;


%%%%量测更新过程
%%Sigma点的更新
I1=sqrt(n+m)*S_next;

 x_sigma_update=x_next;
 for i=2:n+1
     x_sigma_update(:,i)=x_next+I1(:,i-1);
 end
 for i=n+2:L
      x_sigma_update(:,i)=x_next-I1(:,i-n-1);
 end
 %%通过量测方程,得到一系列的粒子
for i=1:n
     y_obser(i)=x_sigma_update(1,i)+x_sigma_update(2,i)*x_sigma_update(3,i);
end
 for i=n+1:L-1
     y_obser(i)=x_sigma_update(1,i-n)+x_sigma_update(2,i-n)*x_sigma_update(3,i-n);
 end
 for i=L;
     y_obser(i)=x_sigma_update(1,i-n-1)+x_sigma_update(2,i-n-1)*x_sigma_update(3,i-n-1);
 end
 %%得到一步预测 
z_next=0;
for i=1:L
     z_next=z_next+w_1*y_obser(i);
end
%%qr分解
for i=1:L
     B(:,i)=sqrt(w)*( y_obser(:,i)-z_next);
end
[R S_zz]=qr([B  sqrt(R)]',0);
%%cholupfdate更新
S_zz=cholupdate(S_zz, y_obser(:,i)-z_next,'-');
 %%得到协方差矩阵Pxz
 Pxz=zeros(3,1);
for i=1:L
   Pxz =Pxz+w*(x_sigma_update(:,i)- x_next)*(y_obser(i)- z_next)';
end
%%%%进行滤波量测更新

%%卡尔曼增益
K = Pxz/(S_zz*S_zz');
U=K*S_zz;
%%估计量的更新
for k=2:81
   x_estimate(:,k)=x_next+K*(z(k-1)-z_next);
end
x_estimate_array= [x_estimate];   %估计值数组

%%方差的更新
S=cholupdate(S_next,U,'-');

 %%误差分析
 for k=1:80
     error(k)=0;
     for j=1:3
     error(k)=error(k)+(x_array(j,k) - x_estimate_array(j,k))^2;
     end
 end
 RMSE=sqrt(error/80);
 %%绘图

 figure
 k=1:1:80;
 plot(k,RMSE,'r-*');xlabel('仿真步长');ylabel('RMSE');

   