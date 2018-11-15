function z_gate= gate_meas_gms(z,gamma,model,m,P)

valid_idx = [];
zlength = size(z,2);
plength = size(m,2);

for j=1:plength
    Sj= model.R + model.H*P(:,:,j)*model.H';
    Vs= chol(Sj); det_Sj= prod(diag(Vs))^2; inv_sqrt_Sj= inv(Vs);
    iSj= inv_sqrt_Sj*inv_sqrt_Sj';
    nu= z- model.H*repmat(m(:,j),[1 zlength]);
    dist= sum((inv_sqrt_Sj'*nu).^2);
    valid_idx= union(valid_idx,find( dist < gamma ));
end
z_gate = z(:,valid_idx);