function pS = compute_pS(model,X)

if isempty(X)
    pS= [];
else
    pS= 0.99*ones(size(X,2),1);
    pS= pS(:);
end
