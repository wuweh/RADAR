%%%=======================================================================
%%% This matlab code implements the square-root version of the CKF (SCKF)
%%% for numerical stability
%%%=======================================================================

clear all;
clc; 
close all;

global Q Qsqrt R Rsqrt;

nExpt = 30;             %No. of Experiments/Trials         
N = 630;                %No. of Time steps            

%%% Process Noise Covariance Q
sigma_processNoise1 = 1e-2;              
sigma_processNoise2 = 1e-1;              
Q = diag([sigma_processNoise1^2  sigma_processNoise2^2 ]); 
Qsqrt = sqrt(Q);        %square-root of Q

%%% Measurement Noise Covariance R
R = 0.005*eye(2);
Rsqrt = sqrt(R);

MSE = zeros(2,N);
estMSE = zeros(2,N);

xestArray = zeros(2,N); 

tic;

[xArray,zArray] = GenerateScenario;

for expt = 1:nExpt
    
    %%%====================================================================
    %%% Initialization
    %%%====================================================================
     
    xkk = [0.3+0.9*rand; pi/2+pi*rand];
    
    Skk = diag([0.9  pi/6]);
    
    fprintf('MC Run in Process = %d\n',expt);      

    for k = 1:N  
        
        %%%================================================================
        %%% By xkk, we mean hat{x}_{k|k}; The corresponding state estimation 
        %%% error covariance Pkk = Skk*Skk'
        %%%================================================================
       %进行qr分解
       [xkk1,Skk1] = Predict(xkk,Skk);
       
       %%%=================================================================
       %%% You may use either Update or Update2; Update2 is a more refined
       %%% version of the square-root CKF
       %%%=================================================================
       
%        [xkk,Skk] = Update(xkk1,Skk1,zArray(:,k));

        [xkk,Skk] = Update2(xkk1,Skk1,zArray(:,k));

        xestArray(:,k) = xkk;

        x = xArray(:,k);
        
        MSE(:,k)= MSE(:,k) + sum((x - xkk).^2); 
        
        estMSE(:,k) = estMSE(:,k)+ trace(Skk*Skk');
        
    end    % time-step
    
end   % expts

toc;

MSE = MSE/(2*nExpt);
estMSE = estMSE/(2*nExpt);

RMSE = MSE.^(0.5);
estRMSE = estMSE.^(0.5);

%%%========================================================================
%%% Plotting
%%%========================================================================

figure;
subplot(2,1,1);
plot(xArray(1,:),'k');
hold on;
plot(xestArray(1,:),'r:');
% xlabel('Time,k','fontsize',16);
ylabel('y(m)','fontsize',16);
legend('Actual','SCKF');
hold off;

subplot(2,1,2);
plot(xArray(2,:),'k');
hold on;
plot(xestArray(2,:),'r:');
xlabel('Time, k','fontsize',16);
ylabel('y(m)','fontsize',16);
legend('Actual','SCKF');
hold off;

figure;
x = 1:630;
% subplot(2,1,1);
semilogy(x,RMSE(1,:),'r');
hold on;
semilogy(x,estRMSE(1,:),'r:');
legend('RMSE','estRMSE');
ylabel('RMSE','fontsize',16);


