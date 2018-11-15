function Z= gen_observation(model,X)
% this observation generation function is for coordinate
% measurements

if isempty(X),
    Z= [];
else
    Z= gen_observation_fn_tgt( model, X, model.D_tgt*randn(size(model.D_tgt,2),size(X,2)) ); %coordinate extraction
end;