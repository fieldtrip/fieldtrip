function [stc] = mne_read_stc_file1(filename)
%
% [stc] = mne_read_stc_file1(filename)
% 
% Reads an stc file. The returned structure has the following fields
%
%     tmin           The first time point of the data in seconds
%     tstep          Time between frames in seconds
%     vertices       vertex indices (1 based)
%     data           The data matrix (nvert * ntime)
%
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
me='MNE:mne_read_stc_file1';
if(nargin ~= 1)
   error(me,'usage: [stc] = mne_read_stc_file1(filename)');
end

try
   stc = mne_read_stc_file(filename);
   stc.vertices = stc.vertices + 1;
catch
   error(me,'%s',mne_omit_first_line(lasterr));
end

return;


