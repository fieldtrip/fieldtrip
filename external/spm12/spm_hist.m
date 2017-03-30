function h = spm_hist(ind,val)
% Generate a weighted histogram - a compiled routine
% FORMAT h = spm_hist(ind,val)
% ind - indices (unsigned byte)
% val - weights
%__________________________________________________________________________
% Copyright (C) 1999-2014 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_hist.m 6157 2014-09-05 18:17:54Z guillaume $


%-This is merely the help file for the compiled routine
%error('spm_hist.c not compiled - see Makefile')

h = accumarray(double(ind(:))+1,double(val(:)),[256 1]);

%-Alternative code for older MATLAB versions
%h = full(sparse(double(ind)+1,ones(size(ind)),val,256,1));
