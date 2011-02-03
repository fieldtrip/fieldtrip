function [s,ds1,ds2,ds3] = spm_sample_priors(b,x1,x2,x3,bg)
% Sample prior probability maps
% FORMAT [s,ds1,ds2,ds3] = spm_sample_priors(b,x1,x2,x3,bg)
% b           - a cell array containing the tissue probability
%               data (see spm_load_priors)
% x1,x2,x3    - coordinates to sample
% bg          - background intensity (i.e. value for points
%               outside FOV)
% s           - sampled values
% ds1,ds2,ds3 - spatial derivatives of sampled values
%____________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$


deg = 3;
lm  = 0;
bg = min(max(bg,lm),(1-lm));
if nargout<=1,
    s      = spm_bsplins(b,x1,x2,x3,[deg deg deg  0 0 0]);
    msk    = find(~isfinite(s));
    s(msk) = bg;
else,
    [s,ds1,ds2,ds3] = spm_bsplins(b,x1,x2,x3,[deg deg deg  0 0 0]);
    msk      = find(~isfinite(s));
    s(msk)   = bg;
    ds1(msk) = 0;
    ds2(msk) = 0;
    ds3(msk) = 0;
end;
