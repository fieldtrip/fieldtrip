function res = nconditions(obj)
% Method for getting the number of unique conditions in the file
% FORMAT res = nconditions(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: nconditions.m 2978 2009-03-27 14:43:08Z guillaume $

try
    res = length(unique(conditions(obj)));
catch
    res = 0;
end