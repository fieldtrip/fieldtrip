function [pth,nam,ext,num] = spm_fileparts(fname)
% Like fileparts, but separates off a comma separated list at the end
% FORMAT [pth,nam,ext,num] = spm_fileparts(fname)
% fname - original filename
% pth   - path
% nam   - filename
% ext   - extension
% num   - comma separated list of values
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$


num = '';
[pth,nam,ext] = fileparts(fname);
ind = find(ext==',');
if ~isempty(ind),
    num = ext(ind(1):end);
    ext = ext(1:(ind(1)-1));
end;

