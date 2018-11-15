function X= gen_newstate_fn(model,Xd,V)

%linear state space equation (CV model)

if ~isnumeric(V)
    if strcmp(V,'noise')
        V= model.sigma_V*model.B*randn(size(model.B,2),size(Xd,2));
    elseif strcmp(V,'noiseless')
        V= zeros(size(model.B,1),size(Xd,2));
    end
end

if isempty(Xd)
    X= [];
else
    X= model.F*Xd+ V;
end