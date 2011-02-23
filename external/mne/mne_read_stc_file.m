function [stc] = mne_read_stc_file(filename)
%
% [stc] = mne_read_stc_file(filename)
% 
% Reads an stc file. The returned structure has the following fields
%
%     tmin           The first time point of the data in seconds
%     tstep          Time between frames in seconds
%     vertices       vertex indices (0 based)
%     data           The data matrix (nvert * ntime)
%
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.7  2006/04/23 15:29:41  msh
%   Added MGH to the copyright
%
%   Revision 1.6  2006/04/10 23:26:54  msh
%   Added fiff reading routines
%
%   Revision 1.5  2005/12/05 20:23:21  msh
%   Added fiff_save_evoked. Improved error handling.
%
%   Revision 1.4  2005/11/21 05:40:24  msh
%   Required number of arguments was wrong.
%
%   Revision 1.3  2005/11/21 03:19:12  msh
%   Improved error handling
%
%   Revision 1.2  2005/11/21 02:15:51  msh
%   Added more routines
%
%   Revision 1.1  2005/11/21 01:41:57  msh
%   Introduced structures and start all function names with mne_
%
%
me='MNE:mne_read_stc_file';
if(nargin ~= 1)
   error(me,'usage: [stc] = mne_read_stc_file(filename)');
end

[fid,message] = fopen(filename,'r','ieee-be');
if fid == -1
    error(me,message)
end

STATUS = fseek(fid, 0, 'eof');
filelen= ftell(fid);
STATUS = fseek(fid, 0, 'bof');

% read tmin in ms
[stc.tmin count]=fread(fid,1,'float32');
stc.tmin=stc.tmin/1000.0;
% read sampling rate in ms
[stc.tstep count]=fread(fid,1,'float32');
stc.tstep=stc.tstep/1000.0;
% read number of vertices/sources
[vertices_n count]=fread(fid,1,'uint32');
% read the source vector
[stc.vertices count]=fread(fid,vertices_n,'uint32');
stc.vertices=uint32(stc.vertices);
% read the number of timepts
[data_n count]=fread(fid,1,'uint32');

if (mod((filelen/4-4-vertices_n),data_n*vertices_n) ~= 0) 
   error(me,'incorrect stc file size')
end


% read the data matrix
[stc.data count]=fread(fid,vertices_n*data_n,'float32');
stc.data=squeeze(reshape(stc.data,vertices_n,data_n));
% close the file
fclose(fid);
return;

