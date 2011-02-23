function fiff_write_id(fid,kind,id)
%
% fiff_write_id(fid,kind,id)
% 
% Writes fiff id 
%
%     fid           An open fif file descriptor
%     kind          The tag kind
%     id            The id to write
%
% If the id argument is missing it will be generated here
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.6  2009/03/30 11:37:37  msh
%   Added copying of measurement info blocks from the original like in mne_browse_raw
%
%   Revision 1.5  2008/06/16 17:27:11  msh
%   Incremented file version number to 1.2
%
%   Revision 1.4  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.3  2006/04/14 11:03:57  msh
%   Changed fiff_write_id write a given id.
%   Added parent id writing.
%
%   Revision 1.2  2006/04/10 23:26:54  msh
%   Added fiff reading routines
%
%   Revision 1.1  2005/12/05 16:01:04  msh
%   Added an initial set of fiff writing routines.
%
%
me='MNE:fiff_write_id';

if nargin == 2
   timezone=5;                              %   Matlab does not know the timezone
   id.version   = bitor(bitshift(1,16),2);  %   Version (1 << 16) | 2
   id.machid(1) = 65536*rand(1);            %   Machine id is random for now
   id.machid(2) = 65536*rand(1);            %   Machine id is random for now
   id.secs      = 3600*(24*(now-datenum(1970,1,1,0,0,0))+timezone);
   id.usecs     = 0;                        %   Do not know how we could get this
elseif nargin ~= 3
   error(me,'Incorrect number of arguments');
end

FIFFT_ID_STRUCT=31;
FIFFV_NEXT_SEQ=0;
%
%
datasize=5*4;                       %   The id comprises five integers
count = fwrite(fid,int32(kind),'int32');
if count ~= 1
    error(me,'write failed');
end
count = fwrite(fid,int32(FIFFT_ID_STRUCT),'int32');
if count ~= 1
    error(me,'write failed');
end
count = fwrite(fid,int32(datasize),'int32');
if count ~= 1
    error(me,'write failed');
end
count = fwrite(fid,int32(FIFFV_NEXT_SEQ),'int32');
if count ~= 1
   error(me,'write failed');
end
%
% Collect the bits together for one write
%
data(1) = id.version;
data(2) = id.machid(1);
data(3) = id.machid(2);
data(4) = id.secs;
data(5) = id.usecs;
count = fwrite(fid,int32(data),'int32');
if count ~= 5
    error(me,'write failed');
end
return;

