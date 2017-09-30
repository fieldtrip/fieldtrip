function trl = ft_trialfun_realtime(cfg)

% FT_TRIALFUN_REALTIME can be used to segment a continuous stream of
% data in real-time. Trials are defined as [begsample endsample offset
% condition]
%
% The configuration structure can contain the following specifications
%   cfg.minsample  = the last sample number that was already considered (passed from rt_process)
%   cfg.blocksize  = in seconds. In case of events, offset is with respect to the trigger.
%   cfg.offset     = the offset wrt the 0 point. In case of no events, offset is wrt
%                    prevSample. E.g., [-0.9 1] will read 1 second blocks with
%                    0.9 second overlap
%   cfg.bufferdata = {'first' 'last'}. If 'last' then only the last block of
%                    interest is read. Otherwise, all well-defined blocks are read (default = 'first')

% Copyright (C) 2009, Marcel van Gerven
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

if ~isfield(cfg,'minsample'),   cfg.minsample = 0;        end
if ~isfield(cfg,'blocksize'),   cfg.blocksize = 0.1;  end
if ~isfield(cfg,'offset'),      cfg.offset = 0; end
if ~isfield(cfg,'bufferdata'),  cfg.bufferdata = 'first'; end
if ~isfield(cfg,'triggers'),    cfg.triggers = [];        end

% blocksize and offset in terms of samples
cfg.blocksize = round(cfg.blocksize * cfg.hdr.Fs);
cfg.offset = round(cfg.offset * cfg.hdr.Fs);

% retrieve trials of interest
if isempty(cfg.event) % asynchronous mode
  trl = trialfun_asynchronous(cfg);
else % synchronous mode
  trl = trialfun_synchronous(cfg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trl = trialfun_asynchronous(cfg)

trl = [];

prevSample = cfg.minsample;

if strcmp(cfg.bufferdata, 'last') % only get last block
  
  % begsample starts blocksize samples before the end
  begsample  = cfg.hdr.nSamples*cfg.hdr.nTrials - cfg.blocksize;
  
  % begsample should be offset samples away from the previous read
  if begsample >= (prevSample + cfg.offset)
    
    endsample  = cfg.hdr.nSamples*cfg.hdr.nTrials;
    
    if begsample < endsample && begsample > 0
      trl = [begsample endsample 0 nan];
    end
  end
  
else % get all blocks
  
  while true
    
    % see whether new samples are available
    newsamples = (cfg.hdr.nSamples*cfg.hdr.nTrials-prevSample);
    
    % if newsamples exceeds the offset plus length specified in blocksize
    if newsamples >= (cfg.offset+cfg.blocksize)
      
      % we do not consider samples < 1
      begsample  = max(1,prevSample+cfg.offset);
      endsample  = max(1,prevSample+cfg.offset+cfg.blocksize);
      
      if begsample < endsample && endsample <= cfg.hdr.nSamples*cfg.hdr.nTrials
        trl = [trl; [begsample endsample 0 nan]];
      end
      prevSample = endsample;
      
    else
      break;
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trl = trialfun_synchronous(cfg)

trl = [];

% process all events
for j=1:length(cfg.event)
  
  if isempty(cfg.triggers)
    curtrig = cfg.event(j).value;
  else
    [m1,curtrig] = ismember(cfg.event(j).value,cfg.triggers);
  end
  
  if isempty(curtrig), curtrig = nan; end
  
  if isempty(cfg.triggers) || (~isempty(m1) && m1)
    % catched a trigger of interest
    
    % we do not consider samples < 1
    begsample = max(1,cfg.event(j).sample + cfg.offset);
    endsample = max(1,begsample + cfg.blocksize);
    
    trl = [trl; [begsample endsample cfg.offset curtrig]];
    
  end
end
