function DB=hbf_BEMOperatorsB_Linear(meshes,coils)
% HBF_BEMOPERATORS_B_LINEAR make DB matrices for linear-basis BEM
%
% DB=HBF_BEMOPERATORS_B_LINEAR(meshes,coils)
%   meshes: BEM geometry, [N x 1] cell array of hbf structs
%   coils:  coil description, hbf struct
%
%   DB: DB matrices in [N x 1] cell array
%
% v160229 (c) Matti Stenroos
Nsurf=length(meshes);
if isfield(coils,'QP')
    fp=coils.QP;
    fpdir=coils.QN;
else
    fp=coils.p;
    fpdir=coils.n;
end
DB=cell(1,Nsurf);
if isfield(coils,'QtoC')
    QtoC=coils.QtoC;
elseif isfield(coils,'QPinds') && ~isempty(coils.QPinds)
    QtoC=QpToCoilsMatrix(coils);
else
    QtoC=1;
end
fprintf('Building %d DB matrices; ',Nsurf);
tic
for I=1:Nsurf
    fprintf('%d ',I);
    DB{I}=QtoC*hbf_DB_Linear(fp,fpdir,meshes{I});  
end
t=toc;
fprintf('done in %.1f seconds.\n',t);

function res=QpToCoilsMatrix(coils)
% function QpToCoils=QpToCoilsMatrix(coils)
coilinds=coils.QPinds;
Nc=size(coilinds,1);
Nqp=coilinds(end,2);
res=zeros(Nc,Nqp);
for I=1:Nc
    inds=(coilinds(I,1):coilinds(I,2));
    w=coils.QW(inds);
    res(I,inds)=w;
end
res=sparse(res);
