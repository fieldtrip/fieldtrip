function [nev] = read_neuralynx_nev(filename, varargin)

% READ_NEURALYNX_NEV reads the event information from the *.nev file in a
% Neuralynx dataset directory
%
% Use as
%   nev = read_neuralynx_hdr(datadir, ...)
%   nev = read_neuralynx_hdr(eventfile, ...)
%
% Optional input arguments should be specified in key-value pairs and may include
%   implementation  should be 1, 2 or 3 (default = 3)
%   value           number or list of numbers
%   mintimestamp    number
%   maxtimestamp    number
%   minnumber       number
%   maxnumber       number
%
% The output structure contains all events and timestamps.

% Copyright (C) 2005, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% get the optional input arguments
flt_value        = ft_getopt(varargin, 'value');
flt_mintimestamp = ft_getopt(varargin, 'mintimestamp');
flt_maxtimestamp = ft_getopt(varargin, 'maxtimestamp');
flt_minnumber    = ft_getopt(varargin, 'minnumber');
flt_maxnumber    = ft_getopt(varargin, 'maxnumber');
implementation   = ft_getopt(varargin, 'implementation', 3);

if ft_filetype(filename, 'neuralynx_ds')
  % replace the directory name by the filename
  if     exist(fullfile(filename, 'Events.Nev'))
    filename = fullfile(filename, 'Events.Nev');
  elseif exist(fullfile(filename, 'Events.nev'))
    filename = fullfile(filename, 'Events.nev');
  elseif exist(fullfile(filename, 'events.Nev'))
    filename = fullfile(filename, 'events.Nev');
  elseif exist(fullfile(filename, 'events.nev'))
    filename = fullfile(filename, 'events.nev');
  end
end

% The file starts with a 16*1024 bytes header in ascii, followed by a
% number of records (c.f. trials).
%
% The format of an event record is
%   int16 PktId
%   int16 PktDataSize
%   int64 TimeStamp
%   int16 EventId
%   int16 TTLValue
%   int16 CRC
%   int32 Dummy
%   int32 Extra[0]
%   ...
%   int32 Extra[7]
%   char EventString[0]
%   ...
%   char EventString[127]
% PktId is usually 0x1002.
% PktDataSize is random data.
% Dummy is random data.
% CRC may contain random data.
% Extra is user-defined data.
% TTLValue is the value sent to the computer on a parallel input port.

fid = fopen_or_error(filename, 'rb', 'ieee-le');

headersize = 16384;
offset     = headersize;
if ~isempty(flt_minnumber)
  offset = offset + (flt_minnumber-1)*184;
end

nev = [];

if implementation==1
  if ~isempty(flt_maxnumber)
    ft_warning('filtering on maximum number not yet implemneted');
  end
  % this is the slow way of reading it
  % it also does not allow filtering
  fseek(fid, offset, 'bof');
  while ~feof(fid)
    nev(end+1).PktStart     = fread(fid, 1, 'int16');
    nev(end  ).PktId        = fread(fid, 1, 'int16');
    nev(end  ).PktDataSize  = fread(fid, 1, 'int16');
    nev(end  ).TimeStamp    = fread(fid, 1, 'uint64');
    nev(end  ).EventId      = fread(fid, 1, 'int16');
    nev(end  ).TTLValue     = fread(fid, 1, 'uint16');
    nev(end  ).CRC          = fread(fid, 1, 'int16');
    nev(end  ).Dummy        = fread(fid, 1, 'int32');
    nev(end  ).Extra        = fread(fid, 8, 'int32');
    nev(end  ).EventString  = fread(fid, 128, 'char');
  end
end

if implementation==2
  if ~isempty(flt_maxnumber)
    ft_warning('filtering on maximum number not yet implemneted');
  end
  % this is a faster way of reading it and it is still using the automatic type conversion from Matlab
  fp = offset;
  fseek(fid, fp+ 0, 'bof'); PktStart       = fread(fid, inf, 'uint16', 184-2);
  num = length(PktStart);
  fseek(fid, fp+ 2, 'bof'); PktId          = fread(fid, num, 'uint16', 184-2);
  fseek(fid, fp+ 4, 'bof'); PktDataSize    = fread(fid, num, 'uint16', 184-2);
  fseek(fid, fp+ 6, 'bof'); TimeStamp      = fread(fid, num, 'uint64=>uint64', 184-8);
  fseek(fid, fp+14, 'bof'); EventId        = fread(fid, num, 'uint16', 184-2);
  fseek(fid, fp+16, 'bof'); TTLValue       = fread(fid, num, 'uint16', 184-2);
  fseek(fid, fp+18, 'bof'); CRC            = fread(fid, num, 'uint16', 184-2);
  fseek(fid, fp+22, 'bof'); Dummy          = fread(fid, num, 'int32' , 184-4);
  % read each of the individual extra int32 values and concatenate them
  fseek(fid, fp+24+0*4, 'bof'); Extra1     = fread(fid, num, 'int32', 184-4);
  fseek(fid, fp+24+1*4, 'bof'); Extra2     = fread(fid, num, 'int32', 184-4);
  fseek(fid, fp+24+2*4, 'bof'); Extra3     = fread(fid, num, 'int32', 184-4);
  fseek(fid, fp+24+3*4, 'bof'); Extra4     = fread(fid, num, 'int32', 184-4);
  fseek(fid, fp+24+4*4, 'bof'); Extra5     = fread(fid, num, 'int32', 184-4);
  fseek(fid, fp+24+5*4, 'bof'); Extra6     = fread(fid, num, 'int32', 184-4);
  fseek(fid, fp+24+6*4, 'bof'); Extra7     = fread(fid, num, 'int32', 184-4);
  fseek(fid, fp+24+7*4, 'bof'); Extra8     = fread(fid, num, 'int32', 184-4);
  Extra = [Extra1 Extra2 Extra3 Extra4 Extra5 Extra6 Extra7 Extra8];
  % read the complete data excluding header as char and cut out the piece with the EventString content
  fseek(fid, fp, 'bof'); EventString        = fread(fid, [184 num], 'char');
  EventString = char(EventString(57:184,:)');
end

if implementation==3
  % this is an even faster way of reading it
  if isempty(flt_minnumber)
    flt_minnumber = 1;
  end
  if isempty(flt_maxnumber)
    flt_maxnumber = inf;
  end
  if fseek(fid, offset, 'bof')~=0
    ft_error(ferror(fid));
  end
  buf = fread(fid, (flt_maxnumber-flt_minnumber+1)*184, 'uint8=>uint8');
  [PktStart, PktId , PktDataSize , TimeStamp , EventId , TTLValue , CRC , Dummy , Extra1 , Extra2 , Extra3 , Extra4 , Extra5 , Extra6 , Extra7 , Extra8 , EventString] = ...
    cstructdecode(buf, 'int16', 'int16', 'int16', 'uint64', 'int16', 'uint16', 'int16', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'char128');
  % the cstructdecode function does not respect the byte order for other-endian numbers
  if bigendian
    PktStart     = swapbytes(PktStart);
    PktId        = swapbytes(PktId);
    PktDataSize  = swapbytes(PktDataSize);
    TimeStamp    = swapbytes(TimeStamp);
    EventId      = swapbytes(EventId);
    TTLValue     = swapbytes(TTLValue);
  % CRC          = swapbytes(CRC);
  % Dummy        = swapbytes(Dummy);
  % Extra1       = swapbytes(Extra1);
  % Extra2       = swapbytes(Extra2);
  % Extra3       = swapbytes(Extra3);
  % Extra4       = swapbytes(Extra4);
  % Extra5       = swapbytes(Extra5);
  % Extra6       = swapbytes(Extra6);
  % Extra7       = swapbytes(Extra7);
  % Extra8       = swapbytes(Extra8);
  end
end

fclose(fid);

if implementation~=1
  % make a selection of events
  sel = true(size(TTLValue));
  if ~isempty(flt_mintimestamp)
    sel = sel & (TimeStamp>=flt_mintimestamp);
  end
  if ~isempty(flt_maxtimestamp)
    sel = sel & (TimeStamp<=flt_maxtimestamp);
  end
  if ~isempty(flt_value)
    tmp = false(size(sel));
    for i=1:length(flt_value)
      tmp = tmp | (TTLValue==flt_value(i));
    end
    sel = sel & tmp;
  end

  % restructure the data into a struct-array, by first making cell-arrays out of it ...
  numsel = sum(sel);
  PktStart      = mat2cell(PktStart(sel)     , ones(1,numsel), 1);
  PktId         = mat2cell(PktId(sel)        , ones(1,numsel), 1);
  PktDataSize   = mat2cell(PktDataSize(sel)  , ones(1,numsel), 1);
  TimeStamp     = mat2cell(TimeStamp(sel)    , ones(1,numsel), 1);
  EventId       = mat2cell(EventId(sel)      , ones(1,numsel), 1);
  TTLValue      = mat2cell(TTLValue(sel)     , ones(1,numsel), 1);
  % CRC           = mat2cell(CRC(sel)          , ones(1,numsel), 1);
  % Dummy         = mat2cell(Dummy(sel)        , ones(1,numsel), 1);
  % Extra1        = mat2cell(Extra1(sel)       , ones(1,numsel), 1);
  % Extra2        = mat2cell(Extra2(sel)       , ones(1,numsel), 1);
  % Extra3        = mat2cell(Extra3(sel)       , ones(1,numsel), 1);
  % Extra4        = mat2cell(Extra4(sel)       , ones(1,numsel), 1);
  % Extra5        = mat2cell(Extra5(sel)       , ones(1,numsel), 1);
  % Extra6        = mat2cell(Extra6(sel)       , ones(1,numsel), 1);
  % Extra7        = mat2cell(Extra7(sel)       , ones(1,numsel), 1);
  % Extra8        = mat2cell(Extra8(sel)       , ones(1,numsel), 1);
  EventString   = mat2cell(EventString(sel,:), ones(1,numsel), 128);
  EventNumber   = mat2cell(find(sel) + flt_minnumber - 1, ones(1,numsel), 1); % this helps in the external filtering for BCI
  % ... and then convert the cell-arrays into a single structure
  nev = struct(...
    'PktStart'    , PktStart     , ...
    'PktId'       , PktId        , ...
    'PktDataSize' , PktDataSize  , ...
    'TimeStamp'   , TimeStamp    , ...
    'EventId'     , EventId      , ...
    'TTLValue'    , TTLValue     , ...
    'EventString' , EventString  , ...
    'EventNumber' , EventNumber);
end

% remove null values and convert to strings
for i=1:length(nev)
  nev(i).EventString = nev(i).EventString(find(nev(i).EventString));
  nev(i).EventString = char(nev(i).EventString(:)');
end

