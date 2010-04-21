function [data, dimord] = data2raw(data);

% DATA2RAW is a helper function that converts various types of averages to
% raw data. This function is used to apply the analysis steps that were
% written for use on preprocessed data also on averaged data.
%
% This function is used in FREQANALYSIS, MEGREALIGN, MEGPLANAR, MEGREPAIR

% Copyright (C) 2005, Robert Oostenveld
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

% determine the type of input data
if isfield(data, 'label') && ~isfield(data, 'avg') && isfield(data, 'trial') && iscell(data.trial)
  dimord = [];                  % output of PREPROCESSING
elseif isfield(data, 'label') && isfield(data, 'avg') && ~isfield(data, 'trial') && ~isfield(data, 'individual')
  dimord = 'chan_time';         % output of TIMELOCKANALYSIS with keeptrials=no or TIMELOCKGRANDAVERAGE with keepindividual=no
elseif isfield(data, 'label') && isfield(data, 'avg') &&  isfield(data, 'trial') && ~isfield(data, 'individual')
  dimord = 'rpt_chan_time';     % output of TIMELOCKANALYSIS with keeptrials=yes
elseif isfield(data, 'label') && isfield(data, 'avg') && ~isfield(data, 'trial') &&  isfield(data, 'individual')
  dimord = 'subj_chan_time';    % output of TIMELOCKGRANDAVERAGE with keepindividual=yes
elseif isfield(data, 'pos') && isfield(data, 'avg') && isfield(data.avg, 'mom')
  % determine the number of source components or orientations
  mom = data.avg.mom{data.inside(1)};
  if size(mom,1)==1
    dimord = 'source_time';
  else
    error('the source activity should be projected on a single orientation');
  end
elseif isfield(data, 'pos') && isfield(data, 'trial') && isfield(data.trial, 'mom')
  % determine the number of source components or orientations
  mom = data.trial(1).mom{data.inside(1)};
  if size(mom,1)==1
    dimord = 'rpt_source_time';
  else
    error('the source activity should be projected on a single orientation');
  end
else
  warning('unrecognized input data');
  dimord = [];
end

if isempty(dimord)
  % nothing to do
  return
end

% convert the input data into raw output, which resembles the output of PREPROCESSING
switch dimord
  case 'chan_time'
    fprintf('converting average to raw data with only a single trial\n');
    data.trial{1} = data.avg;
    data.time     = {data.time};
    data          = rmfield(data, 'avg');

  case 'rpt_chan_time'
    fprintf('converting timelocked trials to raw trials\n');
    tmptrial = {};
    tmptime  = {};
    ntrial = size(data.trial,1);
    nchan  = size(data.trial,2);
    ntime  = size(data.trial,3);
    for i=1:ntrial
      tmptrial{i} = reshape(data.trial(i,:,:), [nchan, ntime]);
      tmptime{i}  = data.time;
    end
    data = rmfield(data, 'trial');
    data.trial = tmptrial;
    data.time  = tmptime;

  case 'subj_chan_time'
    fprintf('converting individual subject averages to raw trials\n');
    tmptrial = {};
    tmptime  = {};
    ntrial = size(data.individual,1);
    nchan  = size(data.individual,2);
    ntime  = size(data.individual,3);
    for i=1:ntrial
      tmptrial{i} = reshape(data.individual(i,:,:), [nchan, ntime]);
      tmptime{i}  = data.time;
    end
    data = rmfield(data, 'individual');
    data.trial = tmptrial;
    data.time  = tmptime;

  case 'source_time'
    fprintf('converting average source reconstructed timecourse to a raw trial\n');
    nsource = length(data.inside);
    ntime = size(data.avg.mom{data.inside(1)},2); % use the first dipole to determine the length
    tmptrial = {};
    tmptime  = {};
    tmptrial{1} = zeros(nsource, ntime);
    for j=1:nsource
      mom = data.avg.mom{data.inside(j)};
      tmptrial{1}(j,:) = mom;
    end
    tmptime{1} = data.time;
    data = rmfield(data, 'avg');
    % add fake raw channel data to the original data structure
    data.trial = tmptrial;
    data.time  = tmptime;
    % add fake channel labels
    data.label = {};
    for i=1:length(data.inside)
      data.label{i} = sprintf('source%d', data.inside(i));
    end
    data.label = data.label(:);

  case 'rpt_source_time'
    fprintf('converting single-trial source reconstructed timecourses to raw trials\n');
    ntrial = length(data.trial);
    nsource = length(data.inside);
    ntime = size(data.trial(1).mom{data.inside(1)},2); % use the first trial and dipole to determine the length
    tmptrial = {};
    tmptime  = {};
    for i=1:ntrial
      tmptrial{i} = zeros(nsource, ntime);
      for j=1:nsource
        mom = data.trial(i).mom{data.inside(j)};
        tmptrial{i}(j,:) = mom;
      end
      tmptime{i} = data.time;
    end
    data = rmfield(data, 'trial');
    % add fake raw channel data to the original data structure
    data.trial = tmptrial;
    data.time  = tmptime;
    % add fake channel labels
    data.label = {};
    for i=1:length(data.inside)
      data.label{i} = sprintf('source%d', data.inside(i));
    end
    data.label = data.label(:);

  otherwise
    % nothing to do

end

% raw data is supposed to contain the sampling frequency
if ~isfield(data, 'fsample')
  % compute the sampling frequency from the first two timepoints
  data.fsample = 1/(data.time{1}(2)-data.time{1}(1));
end

% these fields cannot be used by the external functions that want the data in raw format
if isfield(data, 'avg'),   data = rmfield(data, 'avg');   end
if isfield(data, 'var'),   data = rmfield(data, 'var');   end
if isfield(data, 'cov'),   data = rmfield(data, 'cov');   end
if isfield(data, 'blcov'), data = rmfield(data, 'blcov'); end
