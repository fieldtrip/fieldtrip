function write_plexon_nex(filename, nex)

% WRITE_PLEXON_NEX writes a Plexon *.nex file, which is a file
% containing action-potential (spike) timestamps and waveforms (spike
% channels), event timestamps (event channels), and continuous variable
% data (continuous A/D channels).
%
% Use as
%   write_plexon_nex(filename, nex);
%
% The data structure should contain
%   nex.hdr.FileHeader.Frequency  = TimeStampFreq
%   nex.hdr.VarHeader.Type       = type, 5 for continuous
%   nex.hdr.VarHeader.Name       = label, padded to length 64
%   nex.hdr.VarHeader.WFrequency = sampling rate of continuous channel
%   nex.var.dat                  = data
%   nex.var.ts                   = timestamps
%
% See also READ_PLEXON_NEX, READ_PLEXON_PLX, READ_PLEXON_DDT

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: write_plexon_nex.m,v $
% Revision 1.1  2009/01/14 09:24:45  roboos
% moved even more files from fileio to fileio/privtae, see previous log entry
%
% Revision 1.7  2009/01/08 13:03:09  roboos
% fixed bug in scaling of data, causing it to be clipped (thanks to Conrado)
% fixed bug in seperate scaling of data, the old code was not incorrect, but suboptimal
%
% Revision 1.6  2008/12/15 14:48:22  roboos
% read and write the data in uV/microvolt instead of in mV/milivolt
%
% Revision 1.5  2008/02/04 14:59:00  roboos
% do not rescale the input data if it is int16 already
%
% Revision 1.4  2007/10/08 13:00:24  roboos
% fixed small bug in including cvs identifier as comment in file header
%
% Revision 1.3  2007/03/21 13:01:37  roboos
% reimplemented the calibration and conversion into int16 values
% fixed a bug in the padding of the channel name to 64 chars
%
% Revision 1.2  2007/03/19 17:10:38  roboos
% implemented writing spike waveform data
%
% Revision 1.1  2007/02/21 09:50:51  roboos
% initial implementation
%

% get the optional arguments, these are all required
% FirstTimeStamp = keyval('FirstTimeStamp', varargin);
% TimeStampFreq  = keyval('TimeStampFreq', varargin);

hdr = nex.hdr;

UVtoMV = 1/1000;

switch hdr.VarHeader.Type
  case 5
    dat      = nex.var.dat;                   % this is in microVolt
    buf      = zeros(size(dat), 'int16');
    nchans   = size(dat,1);
    nsamples = size(dat,2);
    nwaves   = 1; % only one continuous datasegment is supported
    if length(hdr.VarHeader)~=nchans
      error('incorrect number of channels');
    end
    % convert the data from floating point into int16 values
    % each channel gets its own optimal calibration factor
    for varlop=1:nchans
      ADMaxValue = double(intmax('int16'));
      ADMaxUV    = max(abs(dat(varlop,:)));   % this is in microVolt
      ADMaxMV    = ADMaxUV/1000;              % this is in miliVolt
      if isa(dat, 'int16')
        % do not rescale data that is already 16 bit
        MVtoAD = 1;
      elseif ADMaxMV==0
        % do not rescale the data if the data is zero
        MVtoAD = 1;
      elseif ADMaxMV>0
        % rescale the data so that it fits into the 16 bits with as little loss as possible
        MVtoAD = ADMaxValue / ADMaxMV;
      end
      buf(varlop,:) = int16(double(dat) * UVtoMV * MVtoAD);
      % remember the calibration value, it should be stored in the variable header
      ADtoMV(varlop) = 1/MVtoAD;
    end
    dat = buf;
    clear buf;

  case 3
    dat      = nex.var.dat;                % this is in microVolt
    nchans   = 1;  % only one channel is supported
    nsamples = size(dat,1);
    nwaves   = size(dat,2);
    if length(hdr.VarHeader)~=nchans
      error('incorrect number of channels');
    end
    % convert the data from floating point into int16 values
    ADMaxValue = double(intmax('int16'));
    ADMaxUV    = max(abs(dat(:)));          % this is in microVolt
    ADMaxMV    = ADMaxUV/1000;              % this is in miliVolt
    if isa(dat, 'int16')
      % do not rescale data that is already 16 bit
      MVtoAD = 1;
    elseif ADMaxMV==0
      % do not rescale the data if the data is zero
      MVtoAD = 1;
    elseif ADMaxMV>0
      % rescale the data so that it fits into the 16 bits with as little loss as possible
      MVtoAD = ADMaxValue / ADMaxMV;
    end
    dat = int16(double(dat) * UVtoMV * MVtoAD);
    % remember the calibration value, it should be stored in the variable header
    ADtoMV = 1/MVtoAD;

  otherwise
    error('unsupported data type')
end % switch type

% determine the first and last timestamp
ts     = nex.var.ts;
ts_beg = min(ts);
ts_end = 0;  % FIXME

fid = fopen(filename, 'wb', 'ieee-le');

% write the file header
write_NexFileHeader;

% write the variable headers
for varlop=1:nchans
  write_NexVarHeader;
end

% write the variable data
for varlop=1:nchans
  write_NexVarData;
end

fclose(fid);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nested function for writing the details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function write_NexFileHeader
    % prepare the two char buffers
    buf1 = padstr('$Id: write_plexon_nex.m,v 1.1 2009/01/14 09:24:45 roboos Exp $', 256);
    buf2 = char(zeros(1, 256));
    % write the stuff to the file
    fwrite(fid, 'NEX1' , 'char');           % NexFileHeader  = string NEX1
    fwrite(fid, 100    , 'int32');          % Version        = version
    fwrite(fid, buf1   , 'char');           % Comment        = comment, 256 bytes
    fwrite(fid, hdr.FileHeader.Frequency, 'double'); % Frequency      = timestamped freq. - tics per second
    fwrite(fid, ts_beg, 'int32');          % Beg            = usually 0, minimum of all the timestamps in the file
    fwrite(fid, ts_end, 'int32');          % End            = maximum timestamp + 1
    fwrite(fid, nchans, 'int32');          % NumVars        = number of variables in the first batch
    fwrite(fid, 0     , 'int32');          % NextFileHeader = position of the next file header in the file, not implemented yet
    fwrite(fid, buf2  , 'char');           % Padding        = future expansion
  end % of the nested function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nested function for writing the details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function write_NexVarHeader
    filheadersize = 544;
    varheadersize = 208;
    offset = filheadersize + nchans*varheadersize + (varlop-1)*nsamples;
    calib = ADtoMV(varlop);
    % prepare the two char buffers
    buf1 = padstr(hdr.VarHeader(varlop).Name, 64);
    buf2 = char(zeros(1, 68));
    % write the continuous variable to the file
    fwrite(fid, hdr.VarHeader.Type, 'int32');        % Type         = 0 - neuron, 1 event, 2- interval, 3 - waveform, 4 - pop. vector, 5 - continuously recorded
    fwrite(fid, 100,                'int32');        % Version      = 100
    fwrite(fid, buf1,               'char');         % Name         = variable name, 1x64 char
    fwrite(fid, offset,             'int32');        % DataOffset   = where the data array for this variable is located in the file
    fwrite(fid, nwaves,             'int32');        % Count        = number of events, intervals, waveforms or weights
    fwrite(fid, 0,                  'int32');        % WireNumber   = neuron only, not used now
    fwrite(fid, 0,                  'int32');        % UnitNumber   = neuron only, not used now
    fwrite(fid, 0,                  'int32');        % Gain         = neuron only, not used now
    fwrite(fid, 0,                  'int32');        % Filter       = neuron only, not used now
    fwrite(fid, 0,                  'double');       % XPos         = neuron only, electrode position in (0,100) range, used in 3D
    fwrite(fid, 0,                  'double');       % YPos         = neuron only, electrode position in (0,100) range, used in 3D
    fwrite(fid, hdr.VarHeader.WFrequency, 'double'); % WFrequency   = waveform and continuous vars only, w/f sampling frequency
    fwrite(fid, calib,     'double');       % ADtoMV       = waveform continuous vars only, coeff. to convert from A/D values to Millivolts
    fwrite(fid, nsamples,  'int32');        % NPointsWave  = waveform only, number of points in each wave
    fwrite(fid, 0,         'int32');        % NMarkers     = how many values are associated with each marker
    fwrite(fid, 0,         'int32');        % MarkerLength = how many characters are in each marker value
    fwrite(fid, buf2,      'char');         % Padding, 1x68 char
  end % of the nested function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nested function for writing the details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function write_NexVarData
    switch hdr.VarHeader.Type
      case 5
        % this code only supports one continuous segment
        index = 0;
        fwrite(fid, ts            , 'int32');       % timestamps, one for each continuous segment
        fwrite(fid, index         , 'int32');       % where to cut the segments, zero offset
        fwrite(fid, dat(varlop,:) , 'int16');       % data
      case 3
        fwrite(fid, ts            , 'int32');       % timestamps, one for each spike
        fwrite(fid, dat           , 'int16');       % waveforms, one for each spike
      otherwise
        error('unsupported data type');
    end % switch
  end % of the nested function

end % of the primary function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction for zero padding a char array to fixed length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = padstr(str, num);
if length(str)>num
  str = str(1:num);
else
  str((end+1):num) = 0;
end
end % of the padstr subfunction
