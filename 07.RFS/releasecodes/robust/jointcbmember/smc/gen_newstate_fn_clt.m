function X= gen_newstate_fn_clt(model,X_old,S,V)
%--- this new state generation function is the random walk model
%but sets the velocity and turn rate components to zero.

if isempty(X_old),
    X= []; return;
end;

X= X_old+ [S;V];



        
   