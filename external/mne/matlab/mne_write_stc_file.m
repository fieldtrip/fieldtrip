function mne_write_stc_file(filename,stc)
%
% mne_write_stc_file(filename,stc)
% 
% writes an stc file
%
%     filename      output file
%     stc           a stucture containing the stc data with fields:
%
%     tmin          The time of the first frame in seconds
%     tstep         Time between frames in seconds
%     vertices      Vertex indices (0 based)
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
%     Revision 1.7  2006/05/05 03:50:40  msh
%     Added routines to compute L2-norm inverse solutions.
%     Added mne_write_inverse_sol_stc to write them in stc files
%     Several bug fixes in other files
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
me='MNE:mne_write_stc_file';
if(nargin ~= 2)
   error(me,'usage: mne_write_stc_file(filename, stc)');
end

[fid,message] = fopen(filename,'w','ieee-be');
if fid == -1
   error(me,message);
end

% write starttime in ms
fwrite(fid,1000*stc.tmin,'float32');
% write sampling rate in ms
fwrite(fid,1000*stc.tstep,'float32');
% write number of vertices
fwrite(fid,length(stc.vertices),'uint32');
% write the vertex indices
for i=1:length(stc.vertices)
    fwrite(fid,stc.vertices(i),'uint32');
end

% write the number of timepts
fwrite(fid,size(stc.data,2),'uint32');
%
% write the data
%
fwrite(fid,stc.data,'float32');

% close the file
fclose(fid);
