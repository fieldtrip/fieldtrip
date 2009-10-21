function [data, dimord] = data2raw(data);

% DATA2RAW is a helper function that converts various types of averages to
% raw data. This function is used to apply the analysis steps that were
% written for use on preprocessed data also on averaged data.
%
% This function is used in FREQANALYSIS, MEGREALIGN, MEGPLANAR, MEGREPAIR

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: data2raw.m,v $
% Revision 1.6  2007/05/16 12:43:31  roboos
% changed a print statement
%
% Revision 1.5  2007/01/17 17:04:17  roboos
% added support for lcmv beamed timecourses, only inside voxels after projectmom
%
% Revision 1.4  2006/02/01 12:26:04  roboos
% made all uses of dimord consistent with the common definition of data dimensions, see the fixdimord() function
%
% Revision 1.3  2005/06/03 07:02:48  roboos
% fixed error for empty dimord in case of raw input data (it should return immediately)
%
% Revision 1.2  2005/06/02 16:00:33  roboos
% removed average from the output structure
%
% Revision 1.1  2005/06/02 12:14:10  roboos
% new implementation for consistent conversion of averaged data (with either keepsubject or keepindividual) to raw trials as they come out of preprocessing
% these two helper functions are from now on used in freqanalysis, megplanar, megrealign, megrepair and in combineplanar
%

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
