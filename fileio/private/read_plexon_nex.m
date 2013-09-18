function [varargout] = read_plexon_nex(filename, varargin)

% READ_PLEXON_NEX reads header or data from a Plexon *.nex file, which
% is a file containing action-potential (spike) timestamps and waveforms
% (spike channels), event timestamps (event channels), and continuous
% variable data (continuous A/D channels).
%
% LFP and spike waveform data that is returned by this function is 
% expressed in microVolt.
%
% Use as
%   [hdr] = read_plexon_nex(filename)
%   [dat] = read_plexon_nex(filename, ...)
%   [dat1, dat2, dat3, hdr] = read_plexon_nex(filename, ...)
%
% Optional arguments should be specified in key-value pairs and can be
%   header      structure with header information
%   feedback    0 or 1
%   tsonly      0 or 1, read only the timestamps and not the waveforms
%   channel     number, or list of numbers (that will result in multiple outputs)
%   begsample   number (for continuous only)
%   endsample   number (for continuous only)
%
% See also READ_PLEXON_PLX, READ_PLEXON_DDT

% Copyright (C) 2007, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

% read header info from file, use Matlabs for automatic byte-ordering
fid = fopen(filename, 'r', 'ieee-le');
fseek(fid, 0, 'eof');
siz = ftell(fid);
fseek(fid, 0, 'bof');

if isempty(hdr)
  if feedback, fprintf('reading header from %s\n', filename); end
  % a NEX file consists of a file header, followed by a number of variable headers
  % sizeof(NexFileHeader) = 544
  % sizeof(NexVarHeader) = 208
  hdr.FileHeader = NexFileHeader(fid);
  if hdr.FileHeader.NumVars<1
    error('no channels present in file');
  end
  hdr.VarHeader = NexVarHeader(fid, hdr.FileHeader.NumVars);
end

for i=1:length(channel)
  chan = channel(i);
  vh = hdr.VarHeader(chan);
  clear buf
  fseek(fid, vh.DataOffset, 'bof');
  switch vh.Type
    case 0
      % Neurons, only timestamps
      buf.ts = fread(fid, [1 vh.Count], 'int32=>int32');

    case 1
      % Events, only timestamps
      buf.ts = fread(fid, [1 vh.Count], 'int32=>int32');

    case 2
      % Interval variables
      buf.begs = fread(fid, [1 vh.Count], 'int32=>int32');
      buf.ends = fread(fid, [1 vh.Count], 'int32=>int32');

    case 3
      % Waveform variables
      buf.ts = fread(fid, [1 vh.Count], 'int32=>int32');
      if ~tsonly
        buf.dat = fread(fid, [vh.NPointsWave vh.Count], 'int16');
        % convert the AD values to miliVolt, subsequently convert from miliVolt to microVolt
        buf.dat = buf.dat * (vh.ADtoMV * 1000);
      end

    case 4
      % Population vector
      error('population vectors are not supported');

    case 5
      % Continuously recorded variables
      buf.ts   = fread(fid, [1 vh.Count], 'int32=>int32');
      buf.indx = fread(fid, [1 vh.Count], 'int32=>int32');
      if vh.Count>1 && (begsample~=1 || endsample~=inf)
        error('reading selected samples from multiple AD segments is not supported');
      end
      if ~tsonly
        numsample = min(endsample - begsample + 1, vh.NPointsWave);
        fseek(fid, (begsample-1)*2, 'cof');
        buf.dat  = fread(fid, [1 numsample], 'int16');
        % convert the AD values to miliVolt, subsequently convert from miliVolt to microVolt
        buf.dat = buf.dat * (vh.ADtoMV * 1000);
      end

    case 6
      % Markers
      buf.ts = fread(fid, [1 vh.Count], 'int32=>int32');
      for j=1:vh.NMarkers
        buf.MarkerNames{j,1} = fread(fid, [1 64], 'uint8=>char');
        for k=1:vh.Count
          buf.MarkerValues{j,k} = fread(fid, [1 vh.MarkerLength], 'uint8=>char');
        end
      end

    otherwise
      error('incorrect channel type');
  end % switch channel type

  % return the data of this channel
  varargout{i} = buf;
end % for channel

% always return the header as last
varargout{end+1} = hdr;

fclose(fid);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hdr = NexFileHeader(fid);
hdr.NexFileHeader  = fread(fid,4,'uint8=>char')';    % string NEX1
hdr.Version        = fread(fid,1,'int32');
hdr.Comment        = fread(fid,256,'uint8=>char')';
hdr.Frequency      = fread(fid,1,'double');         % timestamped freq. - tics per second
hdr.Beg            = fread(fid,1,'int32');          % usually 0
hdr.End            = fread(fid,1,'int32');          % maximum timestamp + 1
hdr.NumVars        = fread(fid,1,'int32');          % number of variables in the first batch
hdr.NextFileHeader = fread(fid,1,'int32');          % position of the next file header in the file, not implemented yet
Padding = fread(fid,256,'uint8=>char')';             % future expansion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hdr = NexVarHeader(fid, numvar);
for varlop=1:numvar
  hdr(varlop).Type         = fread(fid,1,'int32');        % 0 - neuron, 1 event, 2- interval, 3 - waveform, 4 - pop. vector, 5 - continuously recorded
  hdr(varlop).Version      = fread(fid,1,'int32');        % 100
  hdr(varlop).Name         = fread(fid,64,'uint8=>char')'; % variable name
  hdr(varlop).DataOffset   = fread(fid,1,'int32');        % where the data array for this variable is located in the file
  hdr(varlop).Count        = fread(fid,1,'int32');        % number of events, intervals, waveforms or weights
  hdr(varlop).WireNumber   = fread(fid,1,'int32');        % neuron only, not used now
  hdr(varlop).UnitNumber   = fread(fid,1,'int32');        % neuron only, not used now
  hdr(varlop).Gain         = fread(fid,1,'int32');        % neuron only, not used now
  hdr(varlop).Filter       = fread(fid,1,'int32');        % neuron only, not used now
  hdr(varlop).XPos         = fread(fid,1,'double');       % neuron only, electrode position in (0,100) range, used in 3D
  hdr(varlop).YPos         = fread(fid,1,'double');       % neuron only, electrode position in (0,100) range, used in 3D
  hdr(varlop).WFrequency   = fread(fid,1,'double');       % waveform and continuous vars only, w/f sampling frequency
  hdr(varlop).ADtoMV       = fread(fid,1,'double');       % waveform continuous vars only, coeff. to convert from A/D values to Millivolts
  hdr(varlop).NPointsWave  = fread(fid,1,'int32');        % waveform only, number of points in each wave
  hdr(varlop).NMarkers     = fread(fid,1,'int32');        % how many values are associated with each marker
  hdr(varlop).MarkerLength = fread(fid,1,'int32');        % how many characters are in each marker value
  Padding = fread(fid,68,'uint8=>char')';
end
