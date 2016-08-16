function res = delete(this)
% Delete the files of M/EEG dataset from the disk
% FORMAT res = delete(this)
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: delete.m 3378 2009-09-09 16:47:16Z guillaume $

res = 1;

try
    delete(fullfile(path(this), fnamedat(this)));
    delete(fullfile(path(this), fname(this)));
catch
    res = 0;
end
