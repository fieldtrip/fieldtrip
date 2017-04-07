function [data] = raw2data(data, dimord)

% RAW2DATA is a helper function that converts raw data to various types of
% averages. This function is used to apply the analysis steps that were
% written for use on preprocessed data also on averaged data.
%
% This function is the counterpart of DATA2RAW and is used in MEGREALIGN, MEGPLANAR, MEGREPAIR

% Copyright (C) 2005, Robert Oostenveld
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

if isempty(dimord)
  % no conversion is needed
  return;
end

switch dimord
  case 'chan_time'
    fprintf('converting single trial back to average\n');
    data.avg    = data.trial{1};
    data.time   = data.time{1};
    data        = rmfield(data, 'trial');
    data.dimord = dimord;

  case 'rpt_chan_time'
    fprintf('converting raw trials to timelocked trials\n');
    ntrial = length(data.trial);
    nchan  = size(data.trial{1},1);
    ntime  = size(data.trial{1},2);
    tmptrial = zeros(ntrial, nchan, ntime);
    for i=1:ntrial
      tmptrial(i,:,:) = data.trial{i};
    end
    data = rmfield(data, 'trial');
    data.trial  = tmptrial;
    data.avg    = reshape(mean(tmptrial, 1), nchan, ntime);         % recompute the average
    data.var    = reshape(std(tmptrial, [], 1).^2, nchan, ntime);   % recompute the variance
    data.time   = data.time{1};                                     % the time axes of all trials are the same
    data.dimord = dimord;

  case 'subj_chan_time'
    fprintf('converting raw trials to individual subject averages\n');
    nsubj  = length(data.trial);
    nchan  = size(data.trial{1},1);
    ntime  = size(data.trial{1},2);
    tmptrial = zeros(nsubj, nchan, ntime);
    for i=1:nsubj
      tmptrial(i,:,:) = data.trial{i};
    end
    data = rmfield(data, 'trial');
    data.individual = tmptrial;
    data.avg        = reshape(mean(tmptrial, 1), nchan, ntime);         % recompute the average
    data.var        = reshape(std(tmptrial, [], 1).^2, nchan, ntime);   % recompute the variance
    data.time       = data.time{1};                                     % the time axes of all trials are the same
    data.dimord     = dimord;

  otherwise
    warning('unrecognized dimord');
end

