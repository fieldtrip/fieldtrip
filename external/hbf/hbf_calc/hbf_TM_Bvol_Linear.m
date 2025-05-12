function TBvol=hbf_TM_Bvol_Linear(DB,Tphi_full, ci,co)
% HBF_TM_BVOL_LINEAR makes transfer matrix for B_vol
%
%TBvol=HBF_TM_Bvol_LINEAR(DB,Tfull,ci,co)
%   DB:        DB matrix, [N x 1] cell array
%   Tphi_full: BEM transfer matrix for potential
%   ci:        Conductivity inside each BEM surface, [N x 1] 
%   co:        Conductivity outside each BEM surface, [N x 1] 
% 
%   TBvol: Transfer matrix for computing B_vol
%
% v160229 (c) Matti Stenroos

mu0over4pi=1e-7;
Nsurf=length(DB);
startinds=zeros(Nsurf,1);
endinds=zeros(Nsurf,1);
Nrows=size(DB{1},1);
Ncols=0;
for I=1:Nsurf
    startinds(I)=Ncols+1;
    endinds(I)=startinds(I)+size(DB{I},2)-1;
    Ncols=endinds(I);
end

TB=zeros(Nrows,Ncols);
for I=1:Nsurf
          %To make B=Bp+Bvol, make sign different from the convention in some papers and older BEMLIB
          TB(:,startinds(I):endinds(I))=DB{I}*mu0over4pi*(co(I)-ci(I));   
end

TBvol=TB*Tphi_full;
