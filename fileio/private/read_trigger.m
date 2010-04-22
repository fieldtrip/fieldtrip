function [event] = read_trigger(filename, varargin)

% READ_TRIGGER extracts the events from a continuous trigger channel
% This function is a helper function to read_event and can be used for all
% dataformats that have one or multiple continuously sampled TTL channels
% in the data.
%
% The optional trigshift (default is 0) causes the value of the
% trigger to be obtained from a sample that is shifted N samples away
% from the actual flank.
%
% This is a helper function for READ_EVENT
%
% TODO
%  - merge read_ctf_trigger into this function (requires trigshift and bitmasking option)
%  - merge biosemi code into this function (requires bitmasking option)

% Copyright (C) 2008, Robert Oostenveld
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

event = [];

% get the optional input arguments
hdr         = keyval('header',      varargin);
dataformat  = keyval('dataformat',  varargin);
begsample   = keyval('begsample',   varargin);
endsample   = keyval('endsample',   varargin);
chanindx    = keyval('chanindx',    varargin);
detectflank = keyval('detectflank', varargin);
denoise     = keyval('denoise',     varargin); if isempty(denoise),     denoise = 1;      end
trigshift   = keyval('trigshift',   varargin); if isempty(trigshift),   trigshift = 0;    end
trigpadding = keyval('trigpadding', varargin); if isempty(trigpadding), trigpadding = 1;  end
fixctf      = keyval('fixctf',      varargin); if isempty(fixctf),      fixctf = 0;       end
fixneuromag = keyval('fixneuromag', varargin); if isempty(fixneuromag), fixneuromag = 0;  end
fix4dglasgow= keyval('fix4dglasgow', varargin); if isempty(fix4dglasgow), fix4dglasgow = 0; end

if isempty(hdr)
  hdr = ft_read_header(filename);
end

if isempty(begsample)
  begsample = 1;
end

if isempty(endsample)
  endsample = hdr.nSamples*hdr.nTrials;
end

% read the trigger channel as raw data, can safely assume that it is continuous
dat = ft_read_data(filename, 'header', hdr, 'dataformat', dataformat, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', 0);

if isempty(dat)
  % there are no triggers to detect
  return
end

% Detect situations where the channel value changes almost at every time
% step which are likely to be noise
if denoise
  for i=1:length(chanindx)
    if (sum(diff(find(diff(dat(i,:))~=0)) == 1)/length(dat(i,:))) > 0.8
      warning(['trigger channel ' hdr.label{chanindx(i)} ' looks like noise and will be ignored']);
      dat(i,:) = 0;
    end
  end
end

if fixctf
  % correct for reading the data as signed 32-bit integer, whereas it should be interpreted as an unsigned int
  dat(dat<0) = dat(dat<0) + 2^32;
end

if fixneuromag
  % according to Joachim Gross, real events always have triggers > 5
  % this is probably to avoid the noisefloor
  dat(dat<5) = 0;
end

if fix4dglasgow
  % synchronization pulses have a value of 8192 and are set to 0
  dat = dat - bitand(dat, 8192);
  %% triggers containing the first bit assume a value of 4096 when sent by presentation
  %% this does not seem to hold for matlab; check this
  %dat = dat - bitand(dat, 4096)*4095/4096;
end

for i=1:length(chanindx)
  % process each trigger channel independently
  channel = hdr.label{chanindx(i)};
  trig    = dat(i,:);

  if trigpadding
    pad = trig(1);
  else
    pad = 0;
  end

  switch detectflank
    case 'up'
      % convert the trigger into an event with a value at a specific sample
      for j=find(diff([pad trig(:)'])>0)
        event(end+1).type   = channel;
        event(end  ).sample = j + begsample - 1;      % assign the sample at which the trigger has gone down
        event(end  ).value  = trig(j+trigshift);      % assign the trigger value just _after_ going up
      end
    case 'down'
      % convert the trigger into an event with a value at a specific sample
      for j=find(diff([pad trig(:)'])<0)
        event(end+1).type   = channel;
        event(end  ).sample = j + begsample - 1;      % assign the sample at which the trigger has gone down
        event(end  ).value  = trig(j-1-trigshift);    % assign the trigger value just _before_ going down
      end
    case 'both'
      % convert the trigger into an event with a value at a specific sample
      for j=find(diff([pad trig(:)'])>0)
        event(end+1).type   = [channel '_up'];        % distinguish between up and down flank
        event(end  ).sample = j + begsample - 1;      % assign the sample at which the trigger has gone down
        event(end  ).value  = trig(j+trigshift);      % assign the trigger value just _after_ going up
      end
      % convert the trigger into an event with a value at a specific sample
      for j=find(diff([pad trig(:)'])<0)
        event(end+1).type   = [channel '_down'];      % distinguish between up and down flank
        event(end  ).sample = j + begsample - 1;      % assign the sample at which the trigger has gone down
        event(end  ).value  = trig(j-1-trigshift);    % assign the trigger value just _before_ going down
      end
    otherwise
      error('incorrect specification of ''detectflank''');
  end
end
