function l = length(x)
% Overloaded length function for file_array objects
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id: length.m 1143 2008-02-07 19:33:33Z spm $


l = max(size(x));

