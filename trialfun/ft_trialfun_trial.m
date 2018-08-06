function [trl, event] = ft_trialfun_trial(cfg)

% FT_TRIALFUN_TRIAL creates a trial definition that corresponds to the
% events that are returned by FT_READ_EVENT with type='trial'
%
% You can use this function as follows
%   cfg           = [];   
%   cfg.dataset   = string, containing filename or directory
%   cfg.trialfun  = 'ft_trialfun_trial';
%   cfg           = definetrial(cfg);
%   data          = preprocessing(cfg);
%
% See also FT_DEFINETRIAL, FT_PREPROCESSING

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

  tmpsel    = find(smp==smp(sel(i))-offset);
  tmpval    = zeros(1,0);
  for k=1:length(tmpsel)
    if ~isempty(event(tmpsel(k)).value)
      % this assumes per trial that the triggers occur in a standardised
      % order, otherwise the entries per column will have different
      % meanings
      tmpval = [tmpval event(tmpsel(k)).value];
    end
  end
  
  trl(i,1:(length(tmpval)+3))  = [begsample endsample offset tmpval];
end

