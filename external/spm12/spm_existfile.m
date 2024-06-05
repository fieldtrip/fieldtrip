function s = spm_existfile(filename)
% Check if a file exists on disk - a compiled routine
% FORMAT s = spm_existfile(filename)
% filename - filename (can also be a relative or full pathname to a file)
% s        - logical scalar, true if the file exists and false otherwise
%__________________________________________________________________________
%
% This compiled routine is equivalent to:
% >> s = exist(filename,'file') == 2;
% and was written for speed purposes. The differences in behaviour are that
% spm_existfile does not look in MATLAB's search path and does not perform
% tilde '~' expansion.
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


%-This is merely the help file for the compiled routine
%error('spm_existfile.c not compiled - see Makefile')
persistent runOnce
if isempty(runOnce)
    warning('spm_existfile is not compiled for your platform.');
    runOnce = true;
end

s = exist(filename,'file') > 0;
