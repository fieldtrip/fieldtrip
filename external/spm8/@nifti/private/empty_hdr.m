function hdr = empty_hdr
% Create an empty NIFTI-1 header
% FORMAT hdr = empty_hdr
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id: empty_hdr.m 1143 2008-02-07 19:33:33Z spm $


org = niftistruc;
for i=1:length(org),
    hdr.(org(i).label) = org(i).def;
end;

