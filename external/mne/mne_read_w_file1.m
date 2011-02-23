function [w] = mne_read_w_file1(filename)
%
% [w] = mne_read_w_file(filename)
%
% Reads a binary w file into the structure w with the following fields
%
% vertices - vector of vertex indices (1-based)
% data     - vector of data values
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%

me='MNE:mne_read_w_file1';
if(nargin ~= 1)
    error(me,'usage: [w] = mne_read_w_file1(filename)');
end

try
    w = mne_read_w_file(filename);
    w.vertices = w.vertices + 1;
catch
    error(me,'%s',mne_omit_first_line(lasterr));
end

return;
