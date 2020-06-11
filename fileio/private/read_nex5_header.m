function [hdr] = read_nex5_header(filename)

% READ_NEX5_HEADER for Nex Technologies *.nex5 file
%
% Use as
%   [hdr] = read_nex5_header(filename)
%
% See also RAD_NEX5_DATA, READ_NEX5_EVENT

% Copyright (C) 2020, Robert Oostenveld, Alex Kirillov
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

fid = fopen_or_error(filename, 'r', 'ieee-le');

fileheader = Nex5ReadFileHeader(fid);
varheaders = Nex5ReadVarHeaders(fid, fileheader.NumVars);

status = fclose(fid);

% put them together into one struct
hdr.FileHeader = fileheader;
hdr.VarHeader = varheaders;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction to read nex5 file header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hdr = Nex5ReadFileHeader(fid)
hdr.NexFileHeader  = fread(fid,4,'uint8=>char')';   % string NEX5
hdr.Version        = fread(fid,1,'int32');          % valid values are 500, 501 or 502
                                                    % if 500, hdr.End is not specified
                                                    % if 501 or greater, hdr.End is specified
                                                    % if 502 or greater, timestamps can be saved as 64-bit integers
hdr.Comment        = fread(fid,256,'uint8=>char')'; % file comment; UTF-8 encoded
hdr.Frequency      = fread(fid,1,'double');         % timestamped freq. - tics per second; timestamp values are stored in ticks, where tick = 1/Frequency
hdr.Beg            = fread(fid,1,'int64');          % first timestamp of the recording session
hdr.NumVars        = fread(fid,1,'int32');          % number of variables
hdr.MetaStart      = fread(fid,1,'uint64');         % position of metadata at the end of the file
hdr.End            = fread(fid,1,'int64');          % maximum timestamp of the data in the file
Padding = fread(fid,56,'uint8=>char')';             % future expansion
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction to read nex5 variable headers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hdr = Nex5ReadVarHeaders(fid, numvar)
for varindex=1:numvar
  hdr(varindex).Type         = fread(fid,1,'int32');         % 0 - neuron, 1 event, 2- interval, 3 - waveform, 4 - pop. vector, 5 - continuously recorded, 6 - marker
  hdr(varindex).Version      = fread(fid,1,'int32');         % 500
  hdr(varindex).Name         = fread(fid,64,'uint8=>char')'; % variable name; UFT-8 encoded
  hdr(varindex).DataOffset   = fread(fid,1,'uint64');        % where the data array for this variable is located in the file
  hdr(varindex).Count        = fread(fid,1,'uint64');        % number of events, intervals, waveforms or weights
  hdr(varindex).TimestampDataType = fread(fid,1,'int32');    % if 0, timestamps are stored as 32-bit integers; if 1, as 64-bit integers
  hdr(varindex).ContDataType = fread(fid,1,'int32');         % waveforms and continuous variables only; if 0, waveform and continuous values are stored as 16-bit integers; if 1, stored as 32-bit floats
  hdr(varindex).WFrequency   = fread(fid,1,'double');        % waveform and continuous vars only, w/f sampling frequency
  hdr(varindex).Units        = fread(fid,32,'uint8=>char')'; % waveforms and continuous variables only, signal values units, not supported. units are milliVolts
  hdr(varindex).ADtoUnitsCoefficient = fread(fid,1,'double'); % waveforms and continuous variables only, coefficient to convert from A/D values to units.
  hdr(varindex).UnitsOffset = fread(fid,1,'double');         % waveforms and continuous variables only, this offset is used to convert A/D values in units: value_in_units = raw * ADtoUnitsCoefficient + UnitsOffset; ignored if ContinuousDataType == 1
  hdr(varindex).NumberOfDataPoints = fread(fid,1,'uint64');  % waveform variable: number of data points in each wave; continuous variable: overall number of data points in the variable
  hdr(varindex).PrethresholdTimeInSeconds = fread(fid,1,'double'); % waveform variables only, pre-threshold time in seconds
  hdr(varindex).MarkerDataType = fread(fid,1,'int32');       % marker events only; if 0, marker values are stored as strings; if 1, marker values are stored as 32-bit integers
  hdr(varindex).NumberOfMarkerFields = fread(fid,1,'int32'); % marker events only, how many values are associated with each marker
  hdr(varindex).MarkerLength = fread(fid,1,'int32');         % marker events only, how many characters are in each marker value; ignored if MarkerDataType is 1
  hdr(varindex).ContIndexOfFirstPointInFragmentDataType = fread(fid,1,'int32'); % continuous variables only; if 0, indexes are stored as unsigned 32-bit integers; if 1, as unsigned 64-bit integers
  Padding = fread(fid,60,'uint8=>char')';                   % future expansion
end
end

