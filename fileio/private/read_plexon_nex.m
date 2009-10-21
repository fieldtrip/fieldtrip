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
% $Log: read_plexon_nex.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.9  2008/12/15 14:48:22  roboos
% read and write the data in uV/microvolt instead of in mV/milivolt
%
% Revision 1.8  2008/09/30 08:01:04  roboos
% replaced all fread(char=>char) into uint8=>char to ensure that the
% chars are read as 8 bits and not as extended 16 bit characters. The
% 16 bit handling causes problems on some internationalized OS/Matlab
% combinations.
%
% the help of fread specifies "If the precision is 'char' or 'char*1', MATLAB
% reads characters using the encoding scheme associated with the file.
% See FOPEN for more information".
%
% Revision 1.7  2007/10/08 12:59:51  roboos
% give error if no channels present
%
% Revision 1.6  2007/07/19 14:41:56  roboos
% changed indentation and whitespace
%
% Revision 1.5  2007/07/19 08:49:34  roboos
% only give error for multiple continuous segments if specific samples were requested
%
% Revision 1.4  2007/03/26 12:42:20  roboos
% implemented tsonly option to read only the timestamps
% implemented the selection of begin and endsample for continuous channels
%
% Revision 1.3  2007/03/21 12:59:01  roboos
% updated the documentation
% keep timestamps as int32
% convert the AD values to uV for type=3 and 5
% give error instead of warning in case of multiple continuous segments
%
% Revision 1.2  2007/03/14 11:46:16  roboos
% only some whitespace changed
%
% Revision 1.1  2007/01/10 17:28:23  roboos
% new implementation, reusing the code from read_nex_xxx but now with complete support for all known data elements
%

% parse the optional input arguments
hdr       = keyval('header', varargin);
channel   = keyval('channel', varargin);
feedback  = keyval('feedback', varargin);
tsonly    = keyval('tsonly', varargin);
begsample = keyval('begsample', varargin);
endsample = keyval('endsample', varargin);

% set the defaults
if isempty(feedback)
  feedback=0;
end
if isempty(tsonly)
  tsonly=0;
end
if isempty(begsample)
  begsample=1;
end
if isempty(endsample)
  endsample=Inf;
end

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
      ts = fread(fid, [1 vh.Count], 'int32=>int32');
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
