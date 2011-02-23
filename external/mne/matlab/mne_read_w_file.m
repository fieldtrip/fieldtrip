function [w] = mne_read_w_file(filename)
%
% [w] = mne_read_w_file(filename)
%
% Reads a binary w file into the structure w with the following fields
%
% vertices - vector of vertex indices (0-based)
% data     - vector of data values
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%     $Id$
%     
%     Revision 1.6  2006/04/23 15:29:41  msh
%     Added MGH to the copyright
%
%     Revision 1.5  2006/04/10 23:26:54  msh
%     Added fiff reading routines
%
%     Revision 1.4  2005/12/05 20:23:21  msh
%     Added fiff_save_evoked. Improved error handling.
%
%     Revision 1.3  2005/11/21 03:19:12  msh
%     Improved error handling
%
%     Revision 1.2  2005/11/21 02:15:51  msh
%     Added more routines
%
%     Revision 1.1  2005/11/21 01:41:57  msh
%     Introduced structures and start all function names with mne_
%
%
me='MNE:mne_read_w_file';
if(nargin ~= 1)
   error(me,'usage: [w] = mne_read_w_file(filename)');
end

% open it as a big-endian file
[fid,message] = fopen(filename, 'rb', 'b') ;
if (fid < 0)
   error(me,message);
end

fread(fid, 1, 'int16') ;
vnum       = mne_fread3(fid) ;
w.data     = zeros(vnum,1) ;
w.vertices = zeros(vnum,1) ;
for i=1:vnum
  w.vertices(i) = mne_fread3(fid) ; % vertex number (0-based)
  w.data(i)     = fread(fid, 1, 'float') ; % vertex value
end
w.vertices = uint32(w.vertices);

fclose(fid) ;
return;






