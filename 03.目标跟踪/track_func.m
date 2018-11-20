%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   ����������˵����
%   input:  1��N_target��             Ŀ�����
%           2��Target_Info��          Ŀ����Ϣ
%   output: 1��Target_Number��        �ɿ���������
%           2��Target_Info_Output��   �ɿ�������Ϣ
%           3��Target_Number_1��      ��ʱ��������
%           4��Target_Info_Output_1�� ��ʱ������Ϣ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%   Change_History:
%       1��20180609�����ÿ����ֲ��������ޣ������������ڵ�����������������ҹ����㣬�˲�����4ά������׼�������˲�
%       2��20180610���޸Ŀɿ������㼣��ط�����ֱ����XY�᷽�������ж��㼣�Ƿ��뺽�����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [Target_Number, Target_Info_Output,Target_Number_1, Target_Info_Output_1] = track_func(N_target,Target_Info)
    global Track_Table confirm_track_ID temp_track_ID confirm_target_number temp_target_number Conv_P middle_track_ID middle_target_number;    
      
    % ����ȷ�Ϻ���
    if confirm_target_number>0 
        for i=1:confirm_target_number
            deleta_X = abs(Target_Info(1:N_target,1) - Track_Table(confirm_track_ID(i)).X(1));
            deleta_Y = abs(Target_Info(1:N_target,3) - Track_Table(confirm_track_ID(i)).X(3));
            d_2 = (deleta_X.^2 + deleta_Y.^2)/4;
            [value, index] = min(d_2);

            if value<3
                n = 1;
            else
                n = 0;
            end

            if n>0
                 % ɾ���ѹ����ĵ㼣��Ϣ
                Track_Table(confirm_track_ID(i)).M = Target_Info(index,:)';
                Track_Table(confirm_track_ID(i)).FROM = Track_Table(confirm_track_ID(i)).FROM + 1;
                
                %ɾ���ѹ����ĵ㼣��Ϣ
                for k=index:N_target-1 
                    Target_Info(k,:) = Target_Info(k+1,:);
                end         
                N_target = N_target - 1; 
            
                %��ʼ�������˲�
                x = Track_Table(confirm_track_ID(i)).X;
                y = Track_Table(confirm_track_ID(i)).M;
                p = Track_Table(confirm_track_ID(i)).P;               
                [x1,p1] = kalman(x,y,p);                
                Track_Table(confirm_track_ID(i)).X = x1;
                Track_Table(confirm_track_ID(i)).P = p1;
            else 
                 if  Track_Table(confirm_track_ID(i)).no_related == 4 %�����Ĵ�δ������Ŀ��
                    %��ձ�����������Ϣ
                    Track_Table(confirm_track_ID(i)).status = 0;%delete
                    Track_Table(confirm_track_ID(i)).M = zeros(4,1);
                    Track_Table(confirm_track_ID(i)).X = zeros(4,1);
                    Track_Table(confirm_track_ID(i)).Z = zeros(2,1);
                    Track_Table(confirm_track_ID(i)).counter = 0;
                    Track_Table(confirm_track_ID(i)).no_related = 0;
                    Track_Table(confirm_track_ID(i)).Conv = zeros(6,6);
                    Track_Table(confirm_track_ID(i)).FROM = 0;
                    Track_Table(confirm_track_ID(i)).X_AVE = 0;
                else
                    % ����Ԥ��
                    Track_Table(confirm_track_ID(i)).no_related = Track_Table(confirm_track_ID(i)).no_related + 1;
                    Track_Table(confirm_track_ID(i)).M = zeros(4,1);
                    Track_Table(confirm_track_ID(i)).Z(1) = 0;
                    Track_Table(confirm_track_ID(i)).Z(2) = 0; %����δ������Ŀ�꣬��������Ԥ��ֵ��Ϊ��������ֵ��������Ԥ��
                    [Track_Table(confirm_track_ID(i)).X] = linear_filter(Track_Table(confirm_track_ID(i)).X); 
                    Track_Table(confirm_track_ID(i)).FROM = Track_Table(confirm_track_ID(i)).FROM + 1;
                end
            end          
       end
    end
    
       % �����м亽��
     if middle_target_number>0 
        for i=1:middle_target_number
                n = 0;
                deleta_X = abs(Target_Info(1:N_target,1) - Track_Table(middle_track_ID(i)).X(1));
                deleta_Y = abs(Target_Info(1:N_target,3) - Track_Table(middle_track_ID(i)).X(3));
                d_2 = (deleta_X.^2 + deleta_Y.^2)/4;
                [value, index] = min(d_2);
                if value<7
                     n = 1;
                     
                     if value < 3
                       Track_Table(middle_track_ID(i)).kekao_counter = Track_Table(middle_track_ID(i)).kekao_counter + 1;
                     else
                        if Track_Table(middle_track_ID(i)).kekao_counter ~= 0
                            Track_Table(middle_track_ID(i)).kekao_counter =Track_Table(middle_track_ID(i)).kekao_counter - 1;
                        end
                     end
                     
                else
                    n = 0;
                end       
                
            if n==1
                Track_Table(middle_track_ID(i)).no_related = 0;
                Track_Table(middle_track_ID(i)).M = Target_Info(index,:)';
                Track_Table(middle_track_ID(i)).counter = Track_Table(middle_track_ID(i)).counter + 1;
                  for k=index:N_target-1 %ɾ���ѹ����ĵ㼣��Ϣ
                      Target_Info(k,:) = Target_Info(k+1,:);
                  end         
                  N_target = N_target - 1; 
                  
                  if Track_Table(middle_track_ID(i)).kekao_counter == 2  % �������Σ�����Ϊ�ɿ�����
                            Track_Table(middle_track_ID(i)).status = 3;  
                            Track_Table(middle_track_ID(i)).kekao_counter  = 0;
                            Track_Table(middle_track_ID(i)).Middle_counter = 0;
                            Track_Table(middle_track_ID(i)).no_related = 0;
                            
                            Track_Table(middle_track_ID(i)).X = Track_Table(middle_track_ID(i)).M;
                            Track_Table(middle_track_ID(i)).P= Conv_P;
                            Track_Table(middle_track_ID(i)).FROM = 1;
                            %��ʼ�������˲�
                            x = Track_Table(middle_track_ID(i)).X;
                            y = Track_Table(middle_track_ID(i)).M;
                            p = Track_Table(middle_track_ID(i)).P;
                            [Track_Table(middle_track_ID(i)).X, Track_Table(middle_track_ID(i)).P] = kalman(x,y,p);

                  else
                            [Track_Table(middle_track_ID(i)).X] = linear_filter(Track_Table(middle_track_ID(i)).X);
                   end        
            else
                Track_Table(middle_track_ID(i)).no_related  = Track_Table(middle_track_ID(i)).no_related + 1;
                Track_Table(middle_track_ID(i)).counter = 0;
                if Track_Table(middle_track_ID(i)).no_related == 2 %��������û�й����㼣��ɾ���м亽��
                    Track_Table(middle_track_ID(i)).status = 0;%delete
                    Track_Table(middle_track_ID(i)).M = zeros(4,1);
                    Track_Table(middle_track_ID(i)).X = zeros(4,1);
                    Track_Table(middle_track_ID(i)).Z = zeros(2,1);
                    Track_Table(middle_track_ID(i)).Conv = zeros(6,6);
                    Track_Table(middle_track_ID(i)).counter = 0;
                    Track_Table(middle_track_ID(i)).no_related = 0;
                    Track_Table(middle_track_ID(i)).FROM = 0;
                    Track_Table(middle_track_ID(i)).kekao_counter  = 0;
                    Track_Table(middle_track_ID(i)).Middle_counter = 0;
                else
                    if Track_Table(middle_track_ID(i)).kekao_counter ~= 0
                        Track_Table(middle_track_ID(i)).kekao_counter = Track_Table(middle_track_ID(i)).kekao_counter - 1;
                    end
                        % ����Ԥ��
                        [Track_Table(middle_track_ID(i)).X] = linear_filter(Track_Table(middle_track_ID(i)).X);
                end
            end  
        end
     end
     
    % ������ʱ����
     if temp_target_number>0 %������ʱ����
        for i=1:temp_target_number
                deleta_X = abs(Target_Info(1:N_target,1) - Track_Table(temp_track_ID(i)).X(1));
                deleta_Y = abs(Target_Info(1:N_target,3) - Track_Table(temp_track_ID(i)).X(3));
                d_2 = (deleta_X.^2 + deleta_Y.^2)/4;
                [value, index] = min(d_2);
                if value<10
                    n = 1;
                    if value < 3
                       Track_Table(temp_track_ID(i)).kekao_counter = Track_Table(temp_track_ID(i)).kekao_counter + 1;
                    else
                        Track_Table(temp_track_ID(i)).Middle_counter = Track_Table(temp_track_ID(i)).Middle_counter + 1;
                    end
                else
                    n = 0;
                end
                   
            if n==1
                  Track_Table(temp_track_ID(i)).M = Target_Info(index,:)';
                  Track_Table(temp_track_ID(i)).counter = Track_Table(temp_track_ID(i)).counter + 1;
                  for k=index:N_target-1 %ɾ���ѹ����ĵ㼣��Ϣ
                      Target_Info(k,:) = Target_Info(k+1,:);
                  end         
                  N_target = N_target - 1; 
                  if Track_Table(temp_track_ID(i)).counter == 5
                      %�˴��ж�����Ϊ�ɿ����������м亽��
                      if Track_Table(temp_track_ID(i)).kekao_counter >= 3
                            Track_Table(temp_track_ID(i)).status = 3;  %����Ϊ�ɿ�����
                            Track_Table(temp_track_ID(i)).X = Track_Table(temp_track_ID(i)).M;
                            Track_Table(temp_track_ID(i)).P= Conv_P;
                            Track_Table(temp_track_ID(i)).FROM = 1;
                            %��ʼ�������˲�
                            x = Track_Table(temp_track_ID(i)).X;
                            y = Track_Table(temp_track_ID(i)).M;
                            p = Track_Table(temp_track_ID(i)).P;
                            [x1,p1] = kalman(x,y,p);
                            Track_Table(temp_track_ID(i)).X = x1;
                            Track_Table(temp_track_ID(i)).P = p1; 
                            Track_Table(temp_track_ID(i)).kekao_counter  = 0;
                            Track_Table(temp_track_ID(i)).Middle_counter = 0;
                      else
                            Track_Table(temp_track_ID(i)).status = 2;   %����Ϊ�м亽��
                            y = Track_Table(temp_track_ID(i)).M;
                            [x1] = linear_filter(y);
                            Track_Table(temp_track_ID(i)).X = x1;
                            Track_Table(temp_track_ID(i)).kekao_counter  = 0;
                            Track_Table(temp_track_ID(i)).Middle_counter = 0;
                      end        
                  else
                        % ����Ԥ��
                    [Track_Table(temp_track_ID(i)).X] = linear_filter(Track_Table(temp_track_ID(i)).X);
                  end
            end
            if n == 0
                Track_Table(temp_track_ID(i)).status = 0;%delete
                Track_Table(temp_track_ID(i)).M = zeros(4,1);
                Track_Table(temp_track_ID(i)).X = zeros(4,1);
                Track_Table(temp_track_ID(i)).Z = zeros(2,1);
                Track_Table(temp_track_ID(i)).Conv = zeros(6,6);
                Track_Table(temp_track_ID(i)).counter = 0;
                Track_Table(temp_track_ID(i)).no_related = 0;
                Track_Table(temp_track_ID(i)).FROM = 0;
            end  
        end
     end

    % û����ʱ����Ҳû��ȷ�Ϻ���
    if confirm_target_number == 0
        if temp_target_number==0
            if N_target>0
                for i=1:N_target
                    Track_Table(i).status = 1; %1��ʾ��ʱ����
                    Track_Table(i).M = Target_Info(i,:)';
                    Track_Table(i).counter = 1;
                    temp_target_number = temp_target_number + 1;
                    % ����Ԥ��
                    [Track_Table(i).X] = linear_filter(Track_Table(i).M); 
                end
                N_target = 0;
            end
        end
    end
    
    % ��������δ�����ϵĵ�
    if N_target>0
        for i=1:N_target   
            for j=1:64
                if Track_Table(j).status == 0 %��λ�ÿ���
                    Track_Table(j).status = 1;
                    Track_Table(j).M = Target_Info(i,:)';
                    Track_Table(j).counter = 1;
                    % ����Ԥ��
                    [Track_Table(j).X] = linear_filter(Track_Table(j).M); 
                    break;
                end
            end
        end
    end
    N_target = 0;
    
    % ���������
    temp_target_number = 0; 
    temp_track_ID = zeros(256,1);
    middle_target_number = 0;
    middle_track_ID = zeros(256,1);
    confirm_target_number = 0;
    confirm_track_ID = zeros(256,1);
    
    for ii=1:256
        if  Track_Table(ii).status == 1
            temp_target_number = temp_target_number+ 1;
            temp_track_ID(temp_target_number) = ii;       
        end

        if  Track_Table(ii).status == 2
            middle_target_number = middle_target_number+ 1;
            middle_track_ID(middle_target_number) = ii;       
        end

        if  Track_Table(ii).status == 3
            confirm_target_number = confirm_target_number+ 1;
            confirm_track_ID(confirm_target_number) = ii;       
        end
    end
    
    % ������
    Target_Number = confirm_target_number;
    Target_Info_Output = zeros(Target_Number,4);
    for i=1:Target_Number
        Target_Info_Output(i,1:4) = Track_Table(confirm_track_ID(i)).X';
        %��X��Ϣ����ƽ��
        counters = Track_Table(confirm_track_ID(i)).FROM;
        x_value = Track_Table(confirm_track_ID(i)).X_AVE;
        Target_Info_Output(i,1) = (((counters -1)*x_value) + Target_Info_Output(i,1))/counters;
        Track_Table(confirm_track_ID(i)).X_AVE = Target_Info_Output(i,1);
        Target_Info_Output(i,5:8) = Track_Table(confirm_track_ID(i)).M;
    end
    
    Target_Number_1 = temp_target_number;
    Target_Info_Output_1 = zeros(Target_Number_1,4);
    for i=1:Target_Number_1
        Target_Info_Output_1(i,1:4) = Track_Table(temp_track_ID(i)).X';
        Target_Info_Output_1(i,5) = Track_Table(temp_track_ID(i)).counter;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �������˲�����˵����
%   input:  1��X��    ��һ��Ԥ��ֵ
%           2��Z��    ����ֵ
%           3��P��    Э�������
%   output: 1��x��    ��ǰԤ��ֵ
%           2��p��    ����Э�������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function [x,p] = kalman(X,Z,P)
    global Q R ffun hfun;     
    n=numel(X);
    x1 = ffun*X;
    P1 = ffun*P*ffun'+Q;
    z = hfun*x1;
    S = hfun*P1*hfun'+R;
    gain = P1*hfun'/S;
    
    x = x1+gain*(Z-z);
    p = (eye(n)-gain*hfun)*P1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   �����˲�����˵����
%   input:  1��x����ǰ״̬
%   output: 1��X��Ԥ��״̬
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X] = linear_filter(x)  
    global ffun;
    X = ffun*x;
end

            
