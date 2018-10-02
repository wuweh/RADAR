%VD arithmetic
clear;
%%
T=2;
num=50; 
N1=400/T;N2=600/T;N3=610/T;N4=660/T;N5=900/T;
x=zeros(N5,1);
y=zeros(N5,1);
vx=zeros(N5,1);
vy=zeros(N5,1);

%%
x(1)=2000;y(1)=10000;  % initial position
vx(1)=0;vy(1)=-15;     % initial speed
ax=0;ay=0;              %initial acceleration

%%
var=100; %observe noise standard deviation

for i=1:N5-1
    % step 1
    if(i>N1-1&&i<=N2-1)
        ax=0.075;ay=0.075;
        vx(i+1)=vx(i)+ax*T;
        vy(i+1)=vy(i)+ax*T;
     %step 2      
    elseif(i>N2-1&&i<=N3-1)
        ax=0;ay=0;
        vx(i+1)=vx(i);vy(i+1)=vy(i);
    %step 3
    elseif(i>N3-1&&i<=N4-1)
        ax=-0.3;ay=-0.3;
        vx(i+1)=vx(i)+ax*T;
        vy(i+1)=vy(i)+ay*T;
    %step 4    
    else
        ax=0;ay=0;
        vx(i+1)=vx(i);
        vy(i+1)=vy(i);
    end
    %target position info
    x(i+1)=x(i)+vx(i)*T+0.5*ax*T^2;
    y(i+1)=y(i)+vy(i)*T+0.5*ay*T^2;
end

%%



    nx=randn(N5,1); % observe noise
    ny=randn(N5,1); % observe noise
    zx=x+nx;        % observe position
    zy=y+ny;        % observe position
    ki=0;
    low=1;
    high=0;
    u=0;
    ua=0;
    e=0.8;
    xks(1)=zx(1);
    yks(1)=zy(1);
    xks(2)=zx(2);
    yks(2)=zy(2);
    
    o=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1]; %transmit matrix
    h=[1,0,0,0;0,0,1,0]; % observe matrix    
    g=[T^2/2,0;T,0;0,T^2/2;0,T]; %control matrix    
    q=[10000,0;0,10000];  %system covariance   
    %covariance
    perr=[var^2,var^2/T,0,0;var*var/T,2*var^2/(T^2),0,0;0,0,var^2,var^2/T;0,0,var^2/T,2*var^2/(T^2)];
    
    vx=(zx(2)-zx(1))/2;vy=(zy(2)-zy(1))/2;
    xk=[zx(1);vx;zy(1);vy];
    
     for r=3:N5
            z=[zx(r);zy(r)];    % �۲�ֵ
            xk1=o*xk;           % ���µ�ǰʱ�̵�Ԥ��ֵ
            perr1=o*perr*o';    % ���㵱ǰʱ�̵�Э����
            k1=perr1*h'/(h*perr1*h'+q); %���¿���������ϵ��
            xk=xk1+k1*(z-h*xk1);    %���㵱ǰ����ֵ
            perr=(eye(4)-k1*h)*perr1;%���µ�ǰЭ����
            
            %xk���ĸ�Ԫ�طֱ��ǣ�x��vx��y��vy
            xks(r)=xk(1,1);         %��ȡ�������ֵ��x��Ϣ
            yks(r)=xk(3,1);         %��ȡ�������ֵ��y��Ϣ
            vxks(r)=xk(2,1);        %��ȡ�������ֵ��vx��Ϣ
            vyks(r)=xk(4,1);        %��ȡ�������ֵ��vy��Ϣ
            xk1s(r)=xk1(1,1);       %��ȡԤ��ֵ��x��Ϣ
            yk1s(r)=xk1(3,1);       %��ȡԤ��ֵ��y��Ϣ
            vxk1s(r)=xk1(2,1);      %��ȡԤ��ֵ��vx��Ϣ
            vyks1(r)=xk1(4,1);      %��ȡԤ��ֵ��vy��Ϣ
            perr11(r)=perr(1,1);    %��ȡ��ǰЭ�����x��Ϣ
            perr12(r)=perr(1,2);
            perr22(r)=perr(2,2);
     end
     
     figure();
plot(x,y,'k-',zx,zy,'g:',xks,yks,'r-.');
legend('true trace','observation samples','estimated samples');