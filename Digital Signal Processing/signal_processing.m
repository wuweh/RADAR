%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �״��źŴ���Matlab�㷨
% �״������
%       ���δ���200e6��
%       ����90e-6��
%       ����Ƶ�ʣ�2.4096e+10��
%       ����Ƶ�ʣ�1.739MHz��
%       ����������130(��Ч���Զ�60����)��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 20180630��1�����淢��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;
clear all;
close all;
fid = fopen('ADC.txt','r');

C = 3e8;
BW = 200e6;
FC = 2.4096e+10;
lambda = C/FC;
fs = 1.739e6;
T = 90e-6;
kt(2,:) = [2,29,60,128];
dataLength = kt(2,3);
fft_v = 64;

for count = 1:1
    aa=fread(fid,4,'uint32');
    ab=fread(fid,6,'uint8');
    milisecond=fread(fid,1,'uint16');
    ac=fread(fid,2,'uint8');
    ad=fread(fid,4,'uint16');
    
    mouth = ab(1);
    day = ab(2);
    hour = ab(3);
    minute = ab(4);
    second = ab(5);
    
    data_adc1=fread(fid,[aa(4)/2/ac(2),ac(2)],'int16').';
    data_adc2=fread(fid,[aa(4)/2/ac(2),ac(2)],'int16').'; 
end


for m=1:2000
    aa=fread(fid,4,'uint32');
    ab=fread(fid,6,'uint8');
    milisecond=fread(fid,1,'uint16');
    ac=fread(fid,2,'uint8');
    ad=fread(fid,4,'uint16');

    mouth = ab(1);
    day = ab(2);
    hour = ab(3);
    minute = ab(4);
    second = ab(5);
    milisecond = milisecond;

    data_adc1=fread(fid,[aa(4)/2/ac(2),ac(2)],'int16').';
    data_adc2=fread(fid,[aa(4)/2/ac(2),ac(2)],'int16').';
    data1=data_adc1(kt(ac(1)+1,1):end,kt(ac(1)+1,2):(kt(ac(1)+1,2)+kt(ac(1)+1,3)-1));
    data2=data_adc2(kt(ac(1)+1,1):end,kt(ac(1)+1,2):(kt(ac(1)+1,2)+kt(ac(1)+1,3)-1));

    for k=1:fft_v
        adc_avg1=mean(data1(k,:));
        adc_avg2=mean(data2(k,:));
        for k1=1:length(data1(1,:))
            data1(k,k1)=data1(k,k1)-adc_avg1;
            data2(k,k1)=data2(k,k1)-adc_avg2;
        end
    end

    %�ֱ����·ADC��һ��FFT���õ�����ά��������Range_FFT
    Range_FFT_A = fft(data1',kt(ac(1)+1,4))';
    Range_FFT_B = fft(data2',kt(ac(1)+1,4))';

    N_1D_FFT = length(Range_FFT_A(1,:));

    Range_FFT_A=Range_FFT_A.';
    Range_FFT_B=Range_FFT_B.';

    %�Ծ���άFFT����һ��FFT���õ��ٶ�άFFT���
    for i=1:1:N_1D_FFT/2
        fft_2D_data1(i,1:64) = fft([Range_FFT_A(i,:),zeros(1,64-(ac(2)-1))],64);
        show_2D_data1(i,:) = fftshift(fft_2D_data1(i,:)); %ת��Ϊ��0Ƶ��Ϊ�Գ���
        fft_2D_data2(i,1:64) = fft([Range_FFT_B(i,:),zeros(1,64-(ac(2)-1))],64);
        show_2D_data2(i,:) = fftshift(fft_2D_data2(i,:)); %ת��Ϊ��0Ƶ��Ϊ�Գ���
    end

    f2_x_axis = (-64/2:64/2-1);
    AmbData = abs(show_2D_data1);
    AmbData_2 = abs(show_2D_data2); 

    figure(1)
    subplot(2,2,1);
    mesh(f2_x_axis,1:N_1D_FFT/2,abs(AmbData(1:N_1D_FFT/2,1:64)));
    xlabel('V/KHz');
    ylabel('R/KHz');
    zlabel('Amplitude');
    title('ADC1');
   
    subplot(2,2,2);
    mesh(f2_x_axis,1:N_1D_FFT/2,abs(AmbData_2(1:N_1D_FFT/2,1:64)));
    xlabel('V/KHz');
    ylabel('R/KHz');
    zlabel('Amplitude');
    title('ADC2');

    %�ٶ�άCFAR���
    [DotNum,DotLocal] = Speed_CFAR(N_1D_FFT/2,AmbData);

    %����άCFAR���
    AmbData_R = AmbData';
    [targetNum,local1] = Range_CFAR(DotNum, DotLocal,AmbData_R);

    %��ʼ���Ŀ����Ϣ
    row = length(AmbData(:,1));
    col = length(AmbData(1,:));
    result1 = zeros(row,col);
    realTarget1 = 0;
    range = 1;
    speed = 1;

    target_abs = zeros(1,32);
    target_num = 1;
    for n = 1:targetNum
        range = local1((n-1)*2+2);
        speed = local1((n-1)*2+1);
        target_abs(target_num) = AmbData(range,speed);
        target_num = target_num + 1;
    end

    for n=1:targetNum
        range = local1((n-1)*2+2);
        speed = local1((n-1)*2+1);

        if(speed==1)&&((AmbData(range,1)>AmbData(range-1,1))&&(AmbData(range,speed)>AmbData(range+1,speed))&&(AmbData(range,speed)>AmbData(range,speed+1)))
            result1(range,speed) = AmbData(range,speed);
            realLocal1(realTarget1*2+2) = range;
            realLocal1(realTarget1*2+1) = speed;
            realTarget1 = realTarget1 +1;  
        end
        if(speed==64) && ((AmbData(range,64)>AmbData(range-1,64)) && (AmbData(range,speed)>AmbData(range+1,speed)) && (AmbData(range,speed)>AmbData(range,speed-1)))
            result1(range,speed) = AmbData(range,speed);
            realLocal1(realTarget1*2+2) = range;
            realLocal1(realTarget1*2+1) = speed;
            realTarget1 = realTarget1 +1;  
        end
        if(speed>1) && (speed<64) && (AmbData(range,speed)>AmbData(range-1,speed)) && (AmbData(range,speed)>AmbData(range,speed-1)) && (AmbData(range,speed)>AmbData(range+1,speed)) && (AmbData(range,speed)>AmbData(range,speed+1))
            result1(range,speed) = AmbData(range,speed);
            realLocal1(realTarget1*2+2) = range;
            realLocal1(realTarget1*2+1) = speed;
            realTarget1 = realTarget1 +1;  
        end  
    end

    figure(1)
    subplot(2,2,3);
    Y = ((1:N_1D_FFT/2)-1)*fs/kt(2,4)*T*C/BW/2;
    X = (f2_x_axis)*1/T/64*lambda/2; 
    mesh(X,Y,result1(1:N_1D_FFT/2,1:64));
    xlabel('V m/s');
    ylabel('R m');
    zlabel('Amplitude');
    title(realTarget1);

    result1 = zeros(row,col);
    realTarget1 = 0;
    range = 1;
    speed = 1;
    for n = 1:targetNum
        range = local1((n-1)*2+2);
        speed = local1((n-1)*2+1);
        target_abs(target_num) = AmbData(range,speed);
        target_num = target_num + 1;
    end

    for n=1:targetNum
        range = local1((n-1)*2+2);
        speed = local1((n-1)*2+1);
        if(AmbData(range,speed)>AmbData(range-1,speed)) && (AmbData(range,speed)>AmbData(range,speed-1)) && (AmbData(range,speed)>AmbData(range+1,speed))&&(AmbData(range,speed)>AmbData(range,speed+1))
            result1(range,speed) = AmbData(range,speed);
            realLocal1(realTarget1*2+2) = range;
            realLocal1(realTarget1*2+1) = speed;
            realTarget1 = realTarget1 +1;  
        end 
    end

    figure(1)
    subplot(2,2,4);
    Y = ((1:N_1D_FFT/2)-1)*fs/kt(2,4)*T*C/BW/2;
    X = (f2_x_axis)*1/T/64*lambda/2; 
    mesh(X,Y,result1(1:N_1D_FFT/2,1:64));
    xlabel('V m/s');
    ylabel('R m');
    zlabel('Amplitude');
    title(realTarget1);

    %Ŀ����롢�ٶȡ��ǶȽ���
    Range = zeros(32,1);
    Speed = zeros(32,1);
    Azi = zeros(32,1);
    for i=1:realTarget1
        % R = (f_points-1)*fs/FFT_N*T*C/BW/2; 
        Range(i) = (realLocal1((i-1)*2+2)-1)*fs/kt(2,4)*T*C/BW/2;
        % speed = ((1:dopple_fft)-1)*1/T/dopple_fft*lambda/2 (m/s);
        Speed(i) = (realLocal1((i-1)*2+1)-33)*1/T/64*lambda/2;       
        
        %������λ���Ƕ�
        temp_data1 = show_2D_data1(realLocal1((i-1)*2+2),realLocal1((i-1)*2+1))';
        temp_data2 = show_2D_data2(realLocal1((i-1)*2+2),realLocal1((i-1)*2+1));
        temp_data = temp_data1*temp_data2;
        if real(temp_data)>0 && imag(temp_data)>0
            y1 = atan(imag(temp_data)/real(temp_data));
        end
        if real(temp_data)>0 && imag(temp_data)<0
            y1 = atan(imag(temp_data)/real(temp_data));
        end
        if real(temp_data)<0 && imag(temp_data)>0
            y1 = atan(imag(temp_data)/real(temp_data))+pi;
        end
        if real(temp_data)<0 && imag(temp_data)<0
            y1 = atan(imag(temp_data)/real(temp_data))-pi;
        end

        if y1>3.01
            y1 = 3.01;
        end

        if y1<-3.01
            y1 = -3.01;
        end

        y2 = y1*0.01245/2/0.006/pi;
        Azi(i) = asin(y2)/3.14*180+1;
    end
    pause(0.5)
end
fclose(fid);

    
    
    
    
