function D=hbf_BEMOperatorsPhi_LC(meshes)
% HBF_BEMOPERATORS_PHI_LC makes the double-layer (D) matrices for LC BEM
%
% D=HBF_BEMOPERATORS_PHI_LC(meshes)
%   meshes: BEM geometry, [N x 1] cell array of hbf structs
%
%   D: Double-layer matrices (D) in [N x N] cell array
%
% v160229 (c) Matti Stenroos

nsurf=length(meshes);
fprintf('Building %d D matrices; ',nsurf*nsurf);
D=cell(nsurf,nsurf);
indI=repmat(1:nsurf,1,nsurf);
indJ=reshape(indI,nsurf,nsurf)';
indJ=indJ(:)';
tic
for M=1:nsurf*nsurf
         fprintf('%d ',M);
         D{M}=hbf_D_LC(meshes{indI(M)},meshes{indJ(M)});
end
fprintf('done.\n');


for M=1:nsurf
    D{M,M}=SetAutoSolidAngle_LC(D{M,M});
end
toc

function D=SetAutoSolidAngle_LC(D)
% function D=SetAutoSolidAngle_LC(D)
colsums=sum(D,2);
autoangles=-.5-colsums;
for I=1:size(D,1)
    D(I,I)=autoangles(I);
end
