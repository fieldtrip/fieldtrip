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
% $Log: read_trigger.m,v $
% Revision 1.7  2009/05/14 18:53:22  roboos
% bail out immediately if the data is empty
%
% Revision 1.6  2009/02/24 14:25:52  jansch
% added option fix4dglasgow to take the synchronization trigger with value 8192
% out of the trigger-data, prior to flank detection
%
% Revision 1.5  2009/02/09 13:32:36  roboos
% only whitespace
%
% Revision 1.4  2009/01/23 16:18:15  roboos
% changed indentation
%
% Revision 1.3  2009/01/23 12:22:15  vlalit
% A fix to avoid an 'almost infinite' loop in case of noisy event channels.
%
% Revision 1.2  2009/01/23 10:32:55  vlalit
% New reader for Neuromag fif format using the MNE toolbox (http://www.nmr.mgh.harvard.edu/martinos/userInfo/data/sofMNE.php)  implemented by Laurence Hunt.
%
% Revision 1.1  2009/01/14 09:12:16  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.6  2008/05/20 15:12:50  vlalit
% Added trigpadding option to handle channels with baseline different from zero
%
% Revision 1.5  2008/05/15 18:38:53  vlalit
% Fixed the problems with discontinuous files and baseline different than zero
%
% Revision 1.4  2008/05/13 16:48:24  roboos
% added option trigshift (default = 0) for cases where the trigger value should be assigned from a sample not directly after/before the upgoing/downgoing flank
%
% Revision 1.3  2008/05/08 18:32:45  vlalit
% Fixed a bug
%
% Revision 1.2  2008/04/29 14:54:39  roboos
% explicit specification of begsample and endsample, otherwise event.sample remains empty
%
% Revision 1.1  2008/04/29 13:53:50  roboos
% new implementation, works for ctf, bti and neuromag
%

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

if isempty(begsample)
  begsample = 1;
end

if isempty(endsample)
  endsample = hdr.nSamples*hdr.nTrials;
end

% read the trigger channel as raw data, can safely assume that it is continuous
dat = read_data(filename, 'header', hdr, 'dataformat', dataformat, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', 0);

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
