function [s,ds1,ds2,ds3] = spm_sample_priors8(tpm,x1,x2,x3)
% Sample prior probability maps
% FORMAT [s,ds1,ds2,ds3] = spm_sample_priors8(tpm,x1,x2,x3)
% b           - a cell array containing the tissue probability
%               data (see spm_load_priors)
% x1,x2,x3    - coordinates to sample
% s           - sampled values
% ds1,ds2,ds3 - spatial derivatives of sampled values
%
% This function is intended to be used in conjunction with spm_load_priors.
% V = spm_vol(P);
% T = spm_load_priors(V);
% B = spm_sample_priors(T,X,Y,Z);
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_sample_priors8.m 5962 2014-04-17 12:47:43Z spm $

deg  = tpm.deg;
tiny = tpm.tiny;

d  = size(tpm.dat{1});
dx = size(x1);
Kb = numel(tpm.dat);
s  = cell(1,Kb);
msk1 = x1>=1 & x1<=d(1) & x2>=1 & x2<=d(2) & x3>=1 & x3<=d(3);
msk2 = x3<1;
x1 = x1(msk1);
x2 = x2(msk1);
x3 = x3(msk1);
if nargout<=1,
    tot = zeros(dx);
    for k=1:Kb,
        a    = spm_bsplins(tpm.dat{k},x1,x2,x3,[deg deg deg  0 0 0]);
        s{k} = ones(dx)*tpm.bg2(k);
        s{k}(msk1) = exp(a);
        s{k}(msk2) = tpm.bg1(k);
        tot  = tot + s{k};
    end
    msk      = ~isfinite(tot);
    tot(msk) = 1;
    for k=1:Kb,
        s{k}(msk) = tpm.bg2(k);
        s{k}      = s{k}./tot;
    end
else
    ds1 = cell(1,Kb);
    ds2 = cell(1,Kb);
    ds3 = cell(1,Kb);
    tot = zeros(dx);
    for k=1:Kb,
        [a,da1,da2,da3] = spm_bsplins(tpm.dat{k},x1,x2,x3,[deg deg deg  0 0 0]);
        if k==Kb, s{k} = ones(dx); else s{k} = zeros(dx)+tiny; end
        s{k} = ones(dx)*tpm.bg2(k);
        s{k}(msk1) = exp(a);
        s{k}(msk2) = tpm.bg1(k);
        tot    = tot + s{k};
        ds1{k} = zeros(dx); ds1{k}(msk1) = da1;
        ds2{k} = zeros(dx); ds2{k}(msk1) = da2;
        ds3{k} = zeros(dx); ds3{k}(msk1) = da3;
    end
    msk      = find(~isfinite(tot));
    tot(msk) = 1;
    da1      = zeros(dx);
    da2      = zeros(dx);
    da3      = zeros(dx);
    for k=1:Kb,
         s{k}(msk) = tpm.bg1(k);
         s{k}      = s{k}./tot;
         da1       = da1 + s{k}.*ds1{k};
         da2       = da2 + s{k}.*ds2{k};
         da3       = da3 + s{k}.*ds3{k};
    end
    for k=1:Kb,
        ds1{k} = s{k}.*(ds1{k} - da1); ds1{k}(msk) = 0;
        ds2{k} = s{k}.*(ds2{k} - da2); ds2{k}(msk) = 0;
        ds3{k} = s{k}.*(ds3{k} - da3); ds3{k}(msk) = 0;
    end
end
