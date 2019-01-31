function [x_putput, P, likehood] = IMMJPDA_Function(x_pre, P_pre, ffun, hfun, z_noise, num, U)
    global Pg G Q R;

    x_predic_in = ffun*x_pre;
    P_predic = ffun*P_pre*ffun'+G*Q*G'; 
    Z_predic_in = hfun*x_predic_in;
    S = hfun*P_predic*hfun'+ R; 
    K = P_predic*hfun'*inv(S);
    P= P_predic-(1-U(num+1))*K*S*K';

    a=0;         
    b=0;

    x_filter_temp=0;
    for j=1:num 
        x_filter_temp=x_filter_temp+U(j)*(x_predic_in+ K*(z_noise(:,j)- Z_predic_in));
    end
    x_filter_temp=U(j+1)*x_predic_in+x_filter_temp;
    x_filter=x_filter_temp;

    for j=1:num+1
        if j==num+1
            a=x_predic_in;
        else
           a=x_predic_in+ K*(z_noise(:,j)- Z_predic_in);
        end
        b=b+U(j)*(a*a'-x_filter_temp*x_filter_temp');
    end
    P=P+b; 
    x_putput=x_filter;
    
    %计算似然函数
    temp=[];
    for i=1:num
        v1 = z_noise(:,i)-Z_predic_in;
        v = (-0.5*v1'*inv(S)*v1);
        U(i) = exp(v)/(sqrt(2*3.14*det(S)));%似然函数
        LL = U(i)*Pg^(-1);
        temp=[temp LL];
    end
    likehood=sum(temp);
end
