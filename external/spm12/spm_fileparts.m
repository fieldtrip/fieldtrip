function [pth,nam,ext,num] = spm_fileparts(fname)
% Like fileparts, but separates off a comma separated list at the end
% FORMAT [pth,nam,ext,num] = spm_fileparts(fname)
% fname  - original filename
%
% pth    - path
% nam    - filename
% ext    - extension
% num    - comma separated list of values
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2005-2023 Wellcome Centre for Human Neuroimaging


if ~ispc, fname = strrep(fname,'\',filesep); end
[pth,nam,ext] = fileparts(fname);
ind = strfind(ext,',');
if ~isempty(ind)
    if isa(fname,'string') % R2016b onwards
        num = extractAfter(ext,ind(1)-1);
        ext = extractBefore(ext,ind(1));
    else
        num = ext(ind(1):end);
        ext = ext(1:(ind(1)-1));
    end
else
    num = '';
    if isa(fname,'string') % R2016b onwards
        num = string(num);
    end
end
