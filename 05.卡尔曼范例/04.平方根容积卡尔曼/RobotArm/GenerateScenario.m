function [x,y] = GenerateScenario
global Q;

Q1 =Q(1:2,1:2); 
r1 = 0.8; r2 = 0.2;
       noise = Q1 * randn(2,1000);
       theta(1) = 0; theta2(1) = 0; t=1; count = 1;
       for j = 0.3:.1:1.2
           for k = pi/2:.05:1.5*pi
               theta1(t) = j + noise(1,t);
               if  (-1)^count == -1;
                   theta2(t) = k + noise(2,t);  
               else 
                   theta2(t) = 2*pi - k + noise(2,t);
               end    
               xx(t) = r1*cos(theta1(t)) - r2*cos(theta1(t)+theta2(t));
               yy(t) = r1*sin(theta1(t)) - r2*sin(theta1(t)+theta2(t));    
               t = t + 1;
           end
           count = count + 1;
       end   
       x = [theta1;theta2];  % 2-by-630
       y = [xx; yy];   % 2-by-630