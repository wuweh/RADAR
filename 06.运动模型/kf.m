function [x,p] = kf(ffun,X,P,hfun,Z,Q,R)
%KF �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    sizeofdim = size(X,1);          %����ά��
    x1 = ffun*X;                    %����һ��Ԥ�����ֵ
    P1 = ffun*P*ffun'+Q;            %����Ԥ���˲�Э�������
    s = hfun*P1*hfun'+R;
    gain = P1*hfun'/s;  %�����������
    
    x = x1+gain*(Z-hfun*x1);            %�������ֵ
    p =(eye(sizeofdim)-gain*hfun)*P1;   %��Э������·�
end

