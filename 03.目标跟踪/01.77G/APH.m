clc;clear all

A = [ 1 3.8 3 6.6 4.8;
       0.3333 1 0.8667 4.2 3;
       0.3333 1.9333 1 4.6 3;
       0.1543 0.2533 0.2223 1 0.4
       0.2210 0.3333 0.3867 2.8 1;
       ];
   
  for i=1:5
        w1(i)= (A(i,1)*A(i,2)*A(i,3)*A(i,4)*A(i,5))^0.2;
  end
  

w1 = w1';
w2 = w1./(sum(w1))
new =A*w2;
a1 = sum((new(1:5))./(5*w2(1:5)))
 
%  A= [1 3.3 2.8;
%          0.3067 1 0.36;
%          0.4 2.7333 1];
%      
%  for i=1:3
%      w1(i) =  (A(i,1)*A(i,2)*A(i,3))^(1/3);
%  end
% w1 = w1'
% w2 = w1./(sum(w1))
% new =A*w2
% a1 = sum(new(1:3))./(3*w2(1:3))

%  A = [ 1 1/4 1/2;
%             4 1 3;
%             2 1/3 1];
%       for i=1:3
%          w1(i) =  (A(i,1)*A(i,2)*A(i,3))^(1/3);
%      end
% 
%    a1 = w1(1:3)./sum(w1)