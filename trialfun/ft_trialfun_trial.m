function [trl, event] = ft_trialfun_trial(cfg)

% FT_TRIALFUN_TRIAL creates a trial definition that corresponds to the events that
% are returned by FT_READ_EVENT with type='trial'
%
% Use this function by calling
%   [cfg] = ft_definetrial(cfg)
% where the configuration structure should contain
%   cfg.dataset           = string with the filename
%   cfg.trialfun          = 'ft_trialfun_trial'
%
% See also FT_DEFINETRIAL, FT_TRIALFUN_GENERAL

% Copyright (C) 2012, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

if isfield(cfg, 'event')
  % for BCI applications events should be specified in the cfg
  % to prevent reading the same events many times
  event = cfg.event;
else
  event = ft_read_event(cfg.dataset);
end

sel = find(strcmp({event.type}, 'trial'));
trl = zeros(length(sel),3);

smp = [event.sample];

for i=1:length(sel)
  % determine the begin, end and offset for each trial and add it to the Nx3 matrix
  begsample = event(sel(i)).sample;
  endsample = begsample + event(sel(i)).duration - 1;
  offset    = event(sel(i)).offset;

  % this is the value of the trial itself, it might be empty
  value = event(sel(i)).value;

  if isempty(value)
    % try to find corresponding triggers in other events
    extrasel = find(smp==smp(sel(i))-offset);
    value    = nan(1,length(extrasel));
    for k=1:length(extrasel)
      if ~isempty(event(extrasel(k)).value)
        % this assumes per trial that the triggers occur in a standardised
        % order, otherwise the entries per column will have different
        % meanings
        value(k) = event(extrasel(k)).value;
      end
    end
  end

  % store the begin, end and offset, plus the value of the trial
  trl(i,1:(3+length(value))) = [begsample endsample offset value];

end
