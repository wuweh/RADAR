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