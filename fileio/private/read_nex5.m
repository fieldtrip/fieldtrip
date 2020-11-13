function [varargout] = read_nex5(filename, varargin)

% READ_NEX5 reads header or data from a Nex Technologies *.nex5 file, 
% which is a file containing action-potential (spike) timestamps and waveforms
% (spike channels), event timestamps (event channels), and continuous
% variable data (continuous A/D channels).
%
% LFP and spike waveform data that is returned by this function is 
% expressed in microVolt.
%
% Use as
%   [hdr] = read_nex5(filename)
%   [dat] = read_nex5(filename, ...)
%   [dat1, dat2, dat3, hdr] = read_nex5(filename, ...)
%
% Optional arguments should be specified in key-value pairs and can be
%   header      structure with header information
%   feedback    0 or 1
%   tsonly      0 or 1, read only the timestamps and not the waveforms
%   channel     number, or list of numbers (that will result in multiple outputs)
%   begsample   number (for continuous only)
%   endsample   number (for continuous only)
%
% See also READ_NEX5_HEADER
%
% Copyright (C) 2020 Robert Oostenveld, Alex Kirillov
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

% parse the optional input arguments
hdr       = ft_getopt(varargin, 'header');
channel   = ft_getopt(varargin, 'channel');
feedback  = ft_getopt(varargin, 'feedback', false);
tsonly    = ft_getopt(varargin, 'tsonly', false);
begsample = ft_getopt(varargin, 'begsample', 1);
endsample = ft_getopt(varargin, 'endsample', inf);

% start with empty return values and empty data
varargout = {};

if isempty(hdr)
  if feedback, fprintf('reading header from %s\n', filename); end
  hdr = read_nex5_header(filename);
  if hdr.FileHeader.NumVars<1
    ft_error('no channels present in file');
  end
end

% use Matlab for automatic byte-ordering
fid = fopen_or_error(filename, 'r', 'ieee-le');

for i=1:length(channel)
  chan = channel(i);
  vh = hdr.VarHeader(chan);
  clear buf
  status = fseek(fid, vh.DataOffset, 'bof');
  if status < 0;  ft_error('error with fseek');  end
  switch vh.Type
    case 0
      % Neurons, only timestamps
      buf.ts = Nex5ReadTimestamps(fid, vh);
      
    case 1
      % Events, only timestamps
      buf.ts = Nex5ReadTimestamps(fid, vh);

    case 2
      % Interval variables
      buf.begs = Nex5ReadTimestamps(fid, vh);
      buf.ends = Nex5ReadTimestamps(fid, vh);

    case 3
      % Waveform variables
      buf.ts = Nex5ReadTimestamps(fid, vh);
      if ~tsonly
        if vh.ContDataType == 0
          buf.dat = fread(fid, [vh.NumberOfDataPoints vh.Count], 'int16');
        else
          buf.dat = fread(fid, [vh.NumberOfDataPoints vh.Count], 'float32');
        end
        % convert the AD values to miliVolt, subsequently convert from miliVolt to microVolt
        buf.dat = (buf.dat * vh.ADtoUnitsCoefficient + vh.UnitsOffset) * 1000;
      end

    case 4
      % Population vector
      ft_error('population vectors are not supported');

    case 5
      % Continuously recorded variables
      buf.ts = Nex5ReadTimestamps(fid, vh);
      if vh.ContIndexOfFirstPointInFragmentDataType == 0
        buf.indx = fread(fid, [1 vh.Count], 'uint32');
      else
        buf.indx = fread(fid, [1 vh.Count], 'uint64');
      end
      if vh.Count>1 && (begsample~=1 || endsample~=inf)
        ft_error('reading selected samples from multiple AD segments is not supported');
      end
      if ~tsonly
        numsamples = min(endsample - begsample + 1, vh.NumberOfDataPoints);
        status = fseek(fid, (begsample-1)*2, 'cof');
        if status < 0;  ft_error('error with fseek');  end
        if vh.ContDataType == 0
          buf.dat  = fread(fid, [1 numsamples], 'int16');
        else
          buf.dat  = fread(fid, [1 numsamples], 'float32');
        end
        % convert the AD values to miliVolt, subsequently convert from milliVolt to microVolt
        buf.dat = (buf.dat * vh.ADtoUnitsCoefficient + vh.UnitsOffset) * 1000;
      end

    case 6
      % Markers
      buf.ts = Nex5ReadTimestamps(fid, vh);
      for j=1:vh.NumberOfMarkerFields
        buf.MarkerNames{j,1} = fread(fid, [1 64], 'uint8=>char');
        if vh.MarkerDataType == 0
          for k=1:vh.Count
            buf.MarkerValues{j,k} = fread(fid, [1 vh.MarkerLength], 'uint8=>char');
          end
        else
          buf.MarkerValues{j} = fread(fid, vh.Count, 'uint32');
        end
      end

    otherwise
      ft_error('incorrect channel type');
  end % switch channel type

  % return the data of this channel
  varargout{i} = buf;
end % for channel

% always return the header as last
varargout{end+1} = hdr;

fclose(fid);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction to read nex5 timestamps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ts = Nex5ReadTimestamps(fid, varHeader)
if varHeader.TimestampDataType == 0
  ts = fread(fid, [1 varHeader.Count], 'int32');
else
  ts = fread(fid, [1 varHeader.Count], 'int64'); 
end
end
