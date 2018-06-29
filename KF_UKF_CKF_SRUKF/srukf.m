%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THE SQUARE-ROOT UNSCENTED KALMAN FILTER FOR STATE AND PARAMETER-ESTIMATION¡·
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,P]=srukf(ffun,X,P,hfun,Z,Q,R)
    L=numel(X);  
    m=numel(Z);  

    alpha=1e-2;
    ki=0;
    beta=2;

    lambda=alpha^2*(L+ki)-L;            
    c=L+lambda;
    
    Wm=[lambda/c 0.5/c+zeros(1,2*L)]; 
    
    Wc=Wm;
    Wc(1)=Wc(1)+(1-alpha^2+beta);    
    C=sqrt(c);
    
    Xsigmaset=sigmas(X,C,P); 
    
    [X_sigama,Xk_mean,Sk]=ut(ffun,Xsigmaset,Wm,Wc,L,Q); 
    [Y_sigama,Yk_mean,Syk]=ut(hfun,Xsigmaset,Wm,Wc,m,R);   

    Pxy = zeros(6,1);

    for k=1:2*L+1
        Pxy = Pxy + Wc(k)*(X_sigama(:,k)-Xk_mean)*((Y_sigama(:,k)-Yk_mean)');
    end

    kk = Pxy/(Syk')/(Syk);
    X = Xk_mean + kk*(Z-Yk_mean);
    u = kk*Syk;

    for i=1:2
        Sk = cholupdate(Sk,u(:,i),'-');
    end
    P = Sk;
end


function [Xsigma_pre,Mean,Sy]=ut(fun,Xsigma,Wm,Wc,n,COV)
    LL=size(Xsigma,2);
    Mean=zeros(n,1);
    Xsigma_pre=zeros(n,LL);

    for k=1:LL
        Xsigma_pre(:,k)=fun(Xsigma(:,k));
        Mean=Mean+Wm(k)*Xsigma_pre(:,k); 
    end

    A = zeros(n,LL);
    for k=1:2*6+1
        A(:,k)=sqrt(Wc(2))*(Xsigma_pre(:,k)-Mean);
    end

    [~,Sy]=qr([A sqrt(COV)]',0);

    if Wc(1)>0
        Sy=cholupdate(Sy,Xsigma_pre(:,1)-Mean,'+');
    else
        Sy=cholupdate(Sy,Xsigma_pre(:,1)-Mean,'-');
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sigma points around reference point
% Inputs:
% x: reference point
% P: covariance
% c: coefficient
% Output:
% X: Sigma points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xset=sigmas(X,C,P)
    Y = X(:,ones(1,numel(X)));
    Xset = [X Y+C*P Y-C*P];
end