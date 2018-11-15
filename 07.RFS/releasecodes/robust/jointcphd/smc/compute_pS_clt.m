function pS = compute_pS_clt(model,X)

if isempty(X)
    pS= [];
else
    pS= 0.5*ones(size(X,2),1);
end
