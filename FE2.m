function [W] = FE2( vm,dvm,edofMat,K,C,Be,freedofs,p)
ndofs = size(K,1);
nele = size(vm,1);
F = sparse(ndofs,1);
dsp_dvm = sum( vm(:).^p )^(1/p-1);
f = zeros(24,1);
W = zeros(ndofs,1);
for i = 1:nele
    f = -dsp_dvm*( vm(i)^(p-1)*dvm(i,:)*C*Be )';
    F( edofMat(i,:)',1 ) = F( edofMat(i,:)',1 ) + f;
end

 W(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);