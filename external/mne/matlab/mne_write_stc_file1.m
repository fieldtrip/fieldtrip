function mne_write_stc_file1(filename,stc)
%
% mne_write_stc_file1(filename,stc)
% 
% writes an stc file
%
%     filename      output file
%     stc           a stucture containing the stc data with fields:
%
%     tmin          The time of the first frame in seconds
%     tstep         Time between frames in seconds
%     vertices      Vertex indices (1 based)
%     data          The data matrix nvert * ntime
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%     $Id$
%     
me='MNE:mne_write_stc_file1';
if(nargin ~= 2)
   error(me,'usage: mne_write_stc_file1(filename, stc)');
end

stc.vertices = stc.vertices - 1;
try
   mne_write_stc_file(filename,stc);
   stc.vertices = stc.vertices + 1;
catch 
   stc.vertices = stc.vertices - 1; 
   error(me,'%s',mne_omit_first_line(lasterr));
end

return;


