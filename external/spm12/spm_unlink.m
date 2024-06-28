function spm_unlink(varargin)
% Silently delete files on disk - a compiled routine
% FORMAT spm_unlink('file1','file2','file3','file4',...)
%
% Remove the specified file(s) using a system call to unlink().
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 1996-2022 Wellcome Centre for Human Neuroimaging


%-This is merely the help file for the compiled routine
%error('spm_unlink.c not compiled - see Makefile')

rs  = recycle('off');
crs = onCleanup(@() recycle(rs));

ws  = warning('off');
cws = onCleanup(@() warning(ws));

for i=1:numel(varargin)
    delete(varargin{i});
end
