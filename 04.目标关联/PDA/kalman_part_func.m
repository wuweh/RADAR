function [P_predic, Z_predic, S, K] = kalman_part_func(x_in, p_in)
        global C A Q R ;
        P_predic = A*p_in*A'+Q; 
        Z_predic = C*x_in;

        S = C*P_predic*C'+ R; 
        K = P_predic*C'*inv(S);
end