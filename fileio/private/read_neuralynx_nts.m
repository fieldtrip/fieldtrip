function [nts] = read_neuralynx_nts(filename, begrecord, endrecord);

% READ_NEURALYNX_NTS reads spike timestamps
%
% Use as
%   [nts] = read_neuralynx_nts(filename)
%   [nts] = read_neuralynx_nts(filename, begrecord, endrecord)

% Copyright (C) 2006-2007, Robert Oostenveld
%
% $Log: read_neuralynx_nts.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.4  2008/04/29 07:52:31  roboos
% fixed windows related bug
% be consistent with begin and end timestamp in header
%
% Revision 1.3  2008/01/10 12:51:58  roboos
% ensure that it is possible to read only the header, using beg/end = 0/0
%
% Revision 1.2  2007/12/12 16:29:19  roboos
% add first and last timestamp to header, also when no records are read
%
% Revision 1.1  2007/03/21 17:12:04  roboos
% renamed NTE in NTS (filenames and function names)
%
% Revision 1.4  2007/03/21 17:06:57  roboos
% updated the documentation
%
% Revision 1.3  2007/03/21 12:54:37  roboos
% included the 2nd and 3rd input arguments in the function declaration
%
% Revision 1.2  2007/03/19 16:58:08  roboos
% implemented the actual reading of the data from file
%
% Revision 1.1  2006/12/13 15:52:26  roboos
% added empty function, only containing help but no usefull code yet
%

if nargin<2
  begrecord = 1;
end
if nargin<3
  endrecord = inf;
end

% the file starts with a 16*1024 bytes header in ascii, followed by a number of records
hdr = neuralynx_getheader(filename);
fid = fopen(filename, 'rb', 'ieee-le');

% determine the length of the file
fseek(fid, 0, 'eof');
headersize = 16384;
recordsize = 8;
NRecords   = floor((ftell(fid) - headersize)/recordsize);

if begrecord==0 && endrecord==0
  % only read the header  
elseif begrecord<1
  error('cannot read before the first record');
elseif begrecord>NRecords
  error('cannot read beyond the last record')
elseif endrecord>NRecords
  endrecord = NRecords;
end

if NRecords>0
  if (ispc), fclose(fid); end
  % read the timestamp from the first and last record
  hdr.FirstTimeStamp = neuralynx_timestamp(filename, 1);
  hdr.LastTimeStamp  = neuralynx_timestamp(filename, inf);
  if (ispc), fid = fopen(filename, 'rb', 'ieee-le'); end
else
  hdr.FirstTimeStamp = nan;
  hdr.LastTimeStamp  = nan;
end

if begrecord>=1 && endrecord>=begrecord
  % rewind to the first record to be read
  fseek(fid, headersize + (begrecord-1)*recordsize, 'bof');
  numrecord = (endrecord-begrecord+1);
  TimeStamp = fread(fid, numrecord, 'uint64=>uint64');

  nts.TimeStamp = TimeStamp;
end
fclose(fid);

nts.NRecords = NRecords;
nts.hdr      = hdr;
