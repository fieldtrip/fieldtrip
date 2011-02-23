function err = mne_write_w_file(filename, w)
% mne_write_w_file(filename, w)
%
%  writes a binary 'w' file
%
%  filename - name of file to write to
%  w        - a structure with the following fields:
%
% vertices - vector of vertex indices (0-based)
% data     - vector of data values
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
% Revision 1.5  2006/04/23 15:29:41  msh
% Added MGH to the copyright
%
% Revision 1.4  2006/04/10 23:26:54  msh
% Added fiff reading routines
%
% Revision 1.3  2005/12/05 20:23:21  msh
% Added fiff_save_evoked. Improved error handling.
%
% Revision 1.2  2005/11/21 03:19:12  msh
% Improved error handling
%
% Revision 1.1  2005/11/21 02:15:51  msh
% Added more routines
%

me='MNE:mne_write_w_file';
if(nargin ~= 2)
    error(me,'usage: mne_write_w_file(filename, w)');
    return;
end

vnum = length(w.vertices) ;

% open it as a big-endian file
[fid,message] = fopen(filename, 'wb', 'b') ;
if(fid == -1)
    error(me,message);
end

fwrite(fid, 0, 'int16') ;
mne_fwrite3(fid, uint32(vnum)) ;
for i=1:vnum
    mne_fwrite3(fid, w.vertices(i)) ;      % vertex number (0-based)
    fwrite(fid,  w.data(i), 'float') ; % vertex value
end

fclose(fid) ;

return;
