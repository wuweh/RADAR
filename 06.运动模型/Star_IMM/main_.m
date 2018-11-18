clc
clear all
close all
summ=0;
N=10; %总模拟次数

load mytarget 
T=0.1;
R=1;

% %读入Hotarget数据文件
% load Hotarget 
% T=10;
% R=1 %ts xys;
for n=1:N
    y=xys+sqrt(R)*randn(size(xys));
    %%%%%%%%%%%选择模型参数
     a=1/20;  xamax=30;  %当前模型参数
     
     %量测矩阵 start_model
     C_model=[1 0 0];
     C_jerk_model=[1 0 0 0];
    %%%估计横轴
    xe_ca=zeros(3,1);  
    xe_star=zeros(3,1);  
    xe_jerk=zeros(4,1); 
    p=2*eye(3); p_CA = p; p_STAR = p; p_jerk = 2*eye(4);
    %迭代主循环

    x_ca_model=[];    x_star_model=[];    x_jerk_model=[];
    [A_CA,Q_CA]=CAmodel(T,1);

    [A_Jerk,Q_Jerk]=Jerkmodel(T,0.6);
    for i=1:length(y(1,:))
        xa=xe_star(3);
        [A_STAR,Q_STAR,qa]=Starmodel(T,xa,a,xamax);
        [xe_star,p_STAR]=kalmanadfun(A_STAR,C_model,Q_STAR,R,xe_star,y(1,i),p_STAR);
        x_star_model=[x_star_model xe_star];
        
        [xe_jerk,p_jerk]=kalmanadfun(A_Jerk,C_jerk_model,Q_Jerk,R,xe_jerk,y(1,i),p_jerk);
        x_jerk_model=[x_jerk_model xe_jerk];
  
        [xe_ca,p_CA]=kalmanadfun(A_CA,C_model,Q_CA,R,xe_ca,y(1,i),p_CA);
        x_ca_model=[x_ca_model xe_ca]; 
    end
    %%%%估计纵轴
    xe_ca=zeros(3,1);  
    xe_star=zeros(3,1);  
    xe_jerk=zeros(4,1); 
     xe_jerk(1) = y(2,1);
    p=2*eye(3); p_CA = p; p_STAR = p; p_jerk = 2*eye(4);
    y_ca_model=[]; y_star_model = [];y_jerk_model = [];

    [A_CA,Q_CA]=CAmodel(T,1);
    [A_Jerk,Q_Jerk]=Jerkmodel(T,0.6);
    for i=1:length(y(2,:))
        xa=xe_star(3);
        [A_STAR,Q_STAR,qa]=Starmodel(T,xa,a,xamax);
        [xe_star,p_STAR]=kalmanadfun(A_STAR,C_model,Q_STAR,R,xe_star,y(2,i),p_STAR);
        y_star_model=[y_star_model xe_star];
          
        [xe_jerk,p_jerk]=kalmanadfun(A_Jerk,C_jerk_model,Q_Jerk,R,xe_jerk,y(2,i),p_jerk);
        y_jerk_model=[y_jerk_model xe_jerk];

        [xe_ca,p_CA]=kalmanadfun(A_CA,C_model,Q_CA,R,xe_ca,y(2,i),p_CA);
        y_ca_model=[y_ca_model xe_ca];
    end

%     covv=diag(cov(xys'-[C_model*x_ca_model;C_model*y_ca_model]'));
%     summ=summ+covv;
end
% summ/N;
plot(xys(1,:),xys(2,:),'b.-');hold on
plot(C_model*x_ca_model,C_model*y_ca_model,'.-');hold on
plot(C_model*x_star_model,C_model*y_star_model,'.-');hold on
plot(C_jerk_model*x_jerk_model,C_jerk_model*y_jerk_model,'.-');hold off
legend('Real Position','CA Model','Star Model','Jerk Model');

figure 
subplot(2,1,1),
plot(ts,xys(1,:)-C_model*x_star_model);hold on;
plot(ts,xys(1,:)-C_model*x_ca_model);hold on;
plot(ts,xys(1,:)-C_jerk_model*x_jerk_model);hold on;
legend('Star Model','CA Model','Jerk Model');

subplot(2,1,2),
plot(ts,xys(2,:)-C_model*y_star_model);hold on;
plot(ts,xys(2,:)-C_model*y_ca_model);hold on;
plot(ts,xys(2,:)-C_jerk_model*y_jerk_model);hold on;
legend('Star Model','CA Model','Jerk Model');
