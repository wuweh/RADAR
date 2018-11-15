function X= gen_newstate_tg(model,X_old)
%--- this new state generation function follows the linear state spc eqn.
% x= Ax_old + Bv

if isempty(X_old),
    X= [];
else
    mu= X_old(1,:); mu=min(mu,0.999); mu=max(mu,0.001); X_old(1,:)= mu;
    sig= min(0.9*mu.*(1-mu),model.pdvarfac_tg);
    X= gen_newstate_fn_tg(model,X_old,betarnd((mu.*(1-mu)./sig -1).*mu,(mu.*(1-mu)./sig -1).*(1-mu))-mu, model.B*randn(size(model.B,2),size(X_old,2)));
end;