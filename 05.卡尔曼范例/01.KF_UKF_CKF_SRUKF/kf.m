function [x,p] = kf(ffun,X,P,hfun,Z,Q,R)
    n=numel(X);
    x1 = ffun*X;
    P1 = ffun*P*ffun'+Q;
    z = hfun*x1;
    
    S = hfun*P1*hfun'+R;
    gain = P1*hfun'/S;
    
    x = x1+gain*(Z-z);
    p = (eye(n)-gain*hfun)*P1;
end

