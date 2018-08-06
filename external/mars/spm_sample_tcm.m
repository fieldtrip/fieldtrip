function C_slice = spm_sample_tcm(tcm,x1,x2,x3)
% C_slice = spm_sample_tcm(tcm,x1,x2,x3)
%
% Adapted from SPM function spm_sample_priors8(tpm,x1,x2,x3)
% Sample the local version of tissue correlation maps using B-spline
% interpolation.
%
% tcm: a structure containing the TCM data (see spm_load_tcm)
% x1,x2,x3: coordinates to sample
% C_slice: sampled values
%
% This function is intended to be used in conjunction with spm_load_tcm.
% tcm = spm_load_tcm(C);
% C_slice = spm_sample_tcm(tcm,X,Y,Z);
%____________________________________________________________________________
%
% John Ashburner
% Yu (Andy) Huang, 2013-06-04

deg  = tcm.deg;
tiny = eps*eps;

dv = size(tcm.dat);
d  = size(tcm.dat{1});
dx = size(x1);
% Kb = numel(tcm.dat);
C_slice  = cell(dv);

msk1 = x1>=1 & x1<=d(1) & x2>=1 & x2<=d(2) & x3>=1 & x3<=d(3);
msk2 = x3<1;
x1 = x1(msk1);
x2 = x2(msk1);
x3 = x3(msk1);

for k=1:dv(3)
    tot = zeros(dx);
    for j=1:dv(2)
        for i=1:dv(1)
            a = spm_bsplins(tcm.dat{i,j,k},x1,x2,x3,[deg deg deg  0 0 0]); % sampling TCM using B-spline interpolation
                   % deg is NOT same as used by spm_bsplinc???
            C_slice{i,j,k} = ones(dx)*tcm.bg2(i,j,k);
            C_slice{i,j,k}(msk1) = exp(a);
            C_slice{i,j,k}(msk2) = tcm.bg1(i,j,k);
            tot = tot + C_slice{i,j,k};
        end
    end
    msk = ~isfinite(tot);
    tot(msk) = 1;
    for j=1:dv(2)
        for i=1:dv(1)
            C_slice{i,j,k}(msk) = tcm.bg2(i,j,k);
            C_slice{i,j,k} = C_slice{i,j,k}./(tot+tiny);
        end
    end % normalization among the K x K correlation matrix, here is done within one slice of each volume
end