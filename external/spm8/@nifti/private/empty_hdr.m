function hdr = empty_hdr
% Create an empty NIFTI-1 header
% FORMAT hdr = empty_hdr
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id$


org = niftistruc;
for i=1:length(org),
    hdr.(org(i).label) = org(i).def;
end;

