function [data] = raw2data(data, dimord);

% RAW2DATA is a helper function that converts raw data to various types of
% averages. This function is used to apply the analysis steps that were
% written for use on preprocessed data also on averaged data.
%
% This function is the counterpart of DATA2RAW and is used in MEGREALIGN, MEGPLANAR, MEGREPAIR

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: raw2data.m,v $
% Revision 1.5  2006/02/24 15:37:14  roboos
% added the dimord to the output data
%
% Revision 1.4  2006/02/01 12:26:04  roboos
% made all uses of dimord consistent with the common definition of data dimensions, see the fixdimord() function
%
% Revision 1.3  2005/06/30 14:01:52  roboos
% fixed bug in output assignment of avg.time (should not be cell array), thanks to Ingrid
%
% Revision 1.2  2005/06/02 15:29:22  roboos
% fixed error (a variable was called interp instead of data)
%
% Revision 1.1  2005/06/02 12:14:10  roboos
% new implementation for consistent conversion of averaged data (with either keepsubject or keepindividual) to raw trials as they come out of preprocessing
% these two helper functions are from now on used in freqanalysis, megplanar, megrealign, megrepair and in combineplanar
%

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

