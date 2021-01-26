function [event] = read_trigger(filename, varargin)

% READ_TRIGGER extracts the events from a continuous trigger channel
% This function is a helper function to read_event and can be used for all
% dataformats that have one or multiple continuously sampled TTL channels
% in the data.
%
% This is a helper function for FT_READ_EVENT. Please look at the code of
% this function for further details.
%
% TODO
%  - merge read_ctf_trigger into this function (requires trigshift and bitmasking option)
%  - merge biosemi code into this function (requires bitmasking option)
%
% See also FT_READ_EVENT

% Copyright (C) 2008-2020, Robert Oostenveld
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
hdr          = ft_getopt(varargin, 'header'             );
dataformat   = ft_getopt(varargin, 'dataformat'         );
begsample    = ft_getopt(varargin, 'begsample'          );
endsample    = ft_getopt(varargin, 'endsample'          );
chanindx     = ft_getopt(varargin, 'chanindx'           ); % specify -1 in case you don't want to detect triggers
detectflank  = ft_getopt(varargin, 'detectflank'        ); % can be bit, up, down, updiff, downdiff, both
denoise      = ft_getopt(varargin, 'denoise',      true );
trigshift    = ft_getopt(varargin, 'trigshift',    false); % causes the value of the trigger to be obtained from a sample that is shifted N samples away from the actual flank
trigpadding  = ft_getopt(varargin, 'trigpadding',  true );
fixctf       = ft_getopt(varargin, 'fixctf',       false);
fixneuromag  = ft_getopt(varargin, 'fixneuromag',  false);
fix4d8192    = ft_getopt(varargin, 'fix4d8192',    false);
fixbiosemi   = ft_getopt(varargin, 'fixbiosemi',   false);
fixartinis   = ft_getopt(varargin, 'fixartinis',   false);
fixstaircase = ft_getopt(varargin, 'fixstaircase', false);
fixhomer     = ft_getopt(varargin, 'fixhomer',     false);
threshold    = ft_getopt(varargin, 'threshold'          );

if isempty(hdr)
  hdr = ft_read_header(filename);
end

if isempty(begsample)
  begsample = 1;
end

if isempty(endsample)
  endsample = hdr.nSamples*hdr.nTrials;
end

% this is for backward compatibility and can be removed in March 2021
if isequal(detectflank, 'auto')
  % use empty as the default, consistent with how it is done for other options
  detectflank = [];
end

% start with an empty event structure
event = [];

if isempty(chanindx) || isempty(intersect(chanindx, 1:hdr.nChans))
  % there are no triggers to detect
  return
else
  % read the trigger channels as raw data, here we can safely assume that it is continuous
  dat = ft_read_data(filename, 'header', hdr, 'dataformat', dataformat, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', 0);
end

% detect situations where the channel value changes almost at every sample, which are likely to be noise
if istrue(denoise) && isempty(threshold)
  for i=1:length(chanindx)
    % look at how often the value changes, for a clean (i.e. binary) channel this will not be very often
    flanks = find(diff(dat(i,:))~=0);
    if length(flanks) < 0.3 * size(dat,2)
      continue
    end
    % look at the distance between the flanks, for a clean (i.e. binary) channel there will be quite some time between subsequent flanks
    if median(diff(flanks)) > 5
      continue
    end
    % look at the skewness of derivative of the channel, it will be large for a channel with an occasional TTL pulse
    % taking the derivative makes it sensitive for upgoing and downgoing flanks of long TTL pulses
    if skewness(abs(diff(dat(i,:))))>5
      ft_warning(['trigger channel ' hdr.label{chanindx(i)} ' looks like analog TTL pulses and will be thresholded']);
      % use a value halfway the channel-specific extremes
      threshold_value = midrange(dat(i,:));
      dat(i,:) = dat(i,:) >= threshold_value;
    else
      ft_warning(['trigger channel ' hdr.label{chanindx(i)} ' looks like noise and will be ignored']);
      dat(i,:) = 0;
    end
  end
end

if fixbiosemi
  if ft_platform_supports('int32_logical_operations')
    % convert to 32-bit integer representation and only preserve the lowest 24 bits
    dat = bitand(int32(dat), 2^24-1);
    % apparently the 24 bits are still shifted by one byte
    dat = bitshift(dat,-8);
  else
    % find indices of negative numbers
    signbit = find(dat < 0);
    % change type to double (otherwise bitcmp will fail)
    dat = double(dat);
    % make number positive and preserve bits 0-22
    dat(signbit) = bitcmp(abs(dat(signbit))-1,32);
    % apparently the 24 bits are still shifted by one byte
    dat(signbit) = bitshift(dat(signbit),-8);
    % re-insert the sign bit on its original location, i.e. bit24
    dat(signbit) = dat(signbit)+(2^(24-1));
    % typecast the data to ensure that the status channel is represented in 32 bits
    dat = uint32(dat);
  end
  
  byte1 = 2^8  - 1;
  byte2 = 2^16 - 1 - byte1;
  byte3 = 2^24 - 1 - byte1 - byte2;
  
  % get the respective status and trigger bits
  trigger   = bitand(dat, bitor(byte1, byte2)); %  contained in the lower two bytes
  
  % in principle the following bits could also be used, but it would require looking at both flanks for the epoch, cmrange and battery
  % if this code ever needs to be enabled, then it should be done consistently with the biosemi_bdf section in ft_read_event
  % epoch   = int8(bitget(dat, 16+1));
  % cmrange = int8(bitget(dat, 20+1));
  % battery = int8(bitget(dat, 22+1));
  
  % below it will continue with the matrix "dat"
  dat = trigger;
end

if fixctf
  % correct for reading the data as signed 32-bit integer, whereas it should be interpreted as an unsigned int
  dat(dat<0) = dat(dat<0) + 2^32;
end

% fix suggested by Ralph Huonker to deal with triggers that need to be
% interpreted as unsigned integers, rather than signed
if strncmpi(dataformat, 'neuromag', 8) && ~fixneuromag
  for k = 1:size(dat,1)
    switch hdr.chantype{chanindx(1)}
      case 'binary trigger'
        if any(dat(k,:)<0)
          dat(k,:) = double(typecast(int16(dat(k,:)), 'uint16'));
        end
      case 'analog trigger'
        % keep it as it is
      case 'other trigger'
        % keep it as it is
    end
  end
end

if fixneuromag
  % according to Joachim Gross, real events always have triggers > 5
  % this is probably to avoid the noisefloor
  dat(dat<5) = 0;
end

if fix4d8192
  % synchronization pulses have a value of 8192 and are set to 0
  dat = dat - bitand(dat, 8192);
  % triggers containing the first bit assume a value of 4096 when sent by presentation
  % this does not seem to hold for matlab; check this
  % dat = dat - bitand(dat, 4096)*4095/4096;
end

if fixartinis
  % we are dealing with an AD box here, and analog values can be noisy.
  dat = round(10*dat)/10; % steps of 0.1V are to be assumed
end

if fixhomer
  for i=1:numel(chanindx)
    if strcmp(hdr.chantype{chanindx(i)}, 'stimulus')
      % each of the columns of orig.s represents a stimulus type, 1 means on, 0 means off
      % negative values have been editted in Homer and should be ignored
      dat(i,:) = dat(i,:)>0;
    end
  end % for each channel
end

if fixstaircase
  for i=1:numel(chanindx)
    onset  = find(diff([0 dat]>0));
    offset = find(diff([dat 0]<0));
    for j=1:numel(onset)
      % replace all values with the most ocurring value in the window of the TTL pulse
      dat(i,onset:offset) = mode(dat(i,onset:offset));
    end
  end
end

if ~isempty(threshold)
  % the trigger channels contain an analog (and hence noisy) TTL signal and should be thresholded
  if ischar(threshold) % evaluate string (e.g., threshold = 'nanmedian' or 'midrange')
    for i = 1:size(dat,1)
      threshold_value = eval([threshold '(dat(i,:))']);
      % discretize the signal
      dat(i,dat(i,:)< threshold_value) = 0;
      dat(i,dat(i,:)>=threshold_value) = 1;
    end
  else
    % discretize the signal
    dat(dat< threshold) = 0;
    dat(dat>=threshold) = 1;
  end
end

if isempty(dat)
  % either no trigger channels were selected, or no samples
  return
end

if isempty(detectflank)
  % look at the first value in the trigger channel to determine whether the trigger is pulled up or down
  % this fails if the first sample is zero and if the trigger values are negative
  if all(dat(:,1)==0)
    detectflank = 'up';
  else
    detectflank = 'down';
  end
end

for i=1:length(chanindx)
  % process each trigger channel independently
  channel = hdr.label{chanindx(i)};
  trig    = dat(i,:);
  
  if trigpadding
    begpad = trig(1);
    endpad = trig(end);
  else
    begpad = 0;
    endpad = 0;
  end
  
  switch detectflank
    case 'up'
      % convert the trigger into an event with a value at a specific sample
      for j=find(diff([begpad trig])>0)
        event(end+1).type   = channel;
        event(end  ).sample = j + begsample - 1;            % assign the sample at which the trigger has gone up
        event(end  ).value  = trig(j+trigshift);            % assign the trigger value just _after_ going up
      end
    case 'updiff'
      for j=find(diff([begpad trig])>0)
        event(end+1).type   = channel;
        event(end  ).sample = j + begsample - 1;                      % assign the sample at which the trigger has gone up
        event(end  ).value  = trig(j+trigshift)-trig(j+trigshift-1);  % assign the trigger value just _after_ going up minus the value before
      end
    case 'down'
      % convert the trigger into an event with a value at a specific sample
      for j=find(diff([trig endpad])<0)
        event(end+1).type   = channel;
        event(end  ).sample = j + begsample;                % assign the sample at which the trigger has gone down
        event(end  ).value  = trig(j-trigshift);            % assign the trigger value just _before_ going down
      end
    case 'downdiff'
      % convert the trigger into an event with a value at a specific sample
      for j=find(diff([trig endpad])<0)
        event(end+1).type   = channel;
        event(end  ).sample = j + begsample;                          % assign the sample at which the trigger has gone down
        event(end  ).value  = trig(j-trigshift)-trig(j-trigshift+1);  % assign the trigger value just _before_ going up minus the value after
      end
    case 'both'
      % convert the trigger into an event with a value at a specific sample
      difftrace = diff([begpad trig endpad]);
      for j=find(difftrace~=0)
        if difftrace(j)>0
          event(end+1).type   = [channel '_up'];        % distinguish between up and down flank
          event(end  ).sample = j + begsample - 1;      % assign the sample at which the trigger has gone up
          event(end  ).value  = trig(j+trigshift);      % assign the trigger value just _after_ going up
        elseif difftrace(j)<0
          event(end+1).type   = [channel '_down'];      % distinguish between up and down flank
          event(end  ).sample = j + begsample - 1;      % assign the sample at which the trigger has gone down
          event(end  ).value  = trig(j-1-trigshift);    % assign the trigger value just _before_ going down
        end
      end
    case 'any'
      % convert the trigger into an event with a value at a specific sample
      for j=find(diff([begpad trig endpad])~=0)
        event(end+1).type   = channel;
        event(end  ).sample = j + begsample - 1;            % assign the sample at which the trigger has gone up or down
        event(end  ).value  = trig(j+trigshift);            % assign the trigger value just _after_ going up
      end
    case {'bit', 'biton'}
      trig = uint32([begpad trig]);
      for k=1:32
        bitval = bitget(trig, k);                           % get each of the bits separately
        for j=find(~bitval(1:end-1) & bitval(2:end))
          event(end+1).type   = channel;
          event(end  ).sample = j + begsample - 1;          % assign the sample at which the bit has gone up
          event(end  ).value  = 2^(k-1);                    % assign the value represented by this bit
        end % j
      end % k
    case {'bitoff'}
      trig = uint32([trig endpad]);
      for k=1:32
        bitval = bitget(trig, k);                           % get each of the bits separately
        for j=find(bitval(1:end-1) & ~bitval(2:end))
          event(end+1).type   = channel;
          event(end  ).sample = j + begsample;              % assign the sample at which the bit has gone down
          event(end  ).value  = 2^(k-1);                    % assign the value represented by this bit
        end % j
      end % k
    otherwise
      ft_error('incorrect specification of ''detectflank''');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to determine the value that is halfway between the minimum and maximum
% this can be used to threshold, e.g. by specifying this as 'threshold'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = midrange(x)
m = min(x)/2 + max(x)/2;