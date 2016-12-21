function [data] = ft_resampledata(cfg, data)

% FT_RESAMPLEDATA performs a resampling or downsampling of the data
%
% Use as
%   [data] = ft_resampledata(cfg, data)
%
% The data should be organised in a structure as obtained from
% the FT_PREPROCESSING function. The configuration should contain
%   cfg.resamplefs  = frequency at which the data will be resampled (default = 256 Hz)
%   cfg.detrend     = 'no' or 'yes', detrend the data prior to resampling (no default specified, see below)
%   cfg.demean      = 'no' or 'yes', baseline correct the data prior to resampling (default = 'no')
%   cfg.feedback    = 'no', 'text', 'textbar', 'gui' (default = 'text')
%   cfg.trials      = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.sampleindex = 'no' or 'yes', add a channel with the original sample indices (default = 'no')
%
% Instead of specifying cfg.resamplefs, you can also specify a time axis on
% which you want the data to be resampled. This is usefull for merging data
% from two acquisition devides, after resampledata you can call FT_APPENDDATA
% to concatenate the channles from the different acquisition devices.
%   cfg.time        = cell-array with one time axis per trial (i.e. from another dataset)
%   cfg.method      = interpolation method, see INTERP1 (default = 'pchip')
%
% Previously this function used to detrend the data by default. The
% motivation for this is that the data is filtered prior to resampling
% to avoid aliassing and detrending prevents occasional edge artifacts
% of the filters. Detrending is fine for removing slow drifts in data
% priot to frequency analysis, but detrending is not good if you
% subsequenlty want to look at the evoked fields. Therefore the old
% default value 'yes' has been removed. You now explicitely have to
% specify whether you want to detrend (probably so if you want to
% keep your analysis compatible with previous analyses that you did),
% or if you do not want to detrent (recommended in most cases).
% If you observe edge artifacts after detrending, it is recommended
% to apply a baseline correction to the data.
%
% The following fields in the structure 'data' are modified by this function
%   data.fsample
%   data.trial
%   data.time
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_PREPROCESSING

% Copyright (C) 2003-2006, FC Donders Centre, Markus Siegel
% Copyright (C) 2004-2009, FC Donders Centre, Robert Oostenveld
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% ft_checkdata is done further down

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed', {'blc', 'demean'});

% set the defaults
cfg.resamplefs = ft_getopt(cfg, 'resamplefs', []);
cfg.time       = ft_getopt(cfg, 'time',       {});
cfg.detrend    = ft_getopt(cfg, 'detrend',    'no');
cfg.demean     = ft_getopt(cfg, 'demean',     'no');
cfg.feedback   = ft_getopt(cfg, 'feedback',   'text');
cfg.trials     = ft_getopt(cfg, 'trials',     'all', 1);
cfg.method     = ft_getopt(cfg, 'method',     'pchip');
cfg.sampleindex = ft_getopt(cfg, 'sampleindex', 'no');

% give the user control over whether to use resample (applies anti-aliasing
% filter) or downsample (does not apply filter)
cfg.resamplemethod = ft_getopt(cfg, 'resamplemethod', 'resample');

% store original datatype
convert = ft_datatype(data);

% check if the input data is valid for this function, this will convert it to raw if needed
data = ft_checkdata(data, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes');

%set default resampling frequency
if isempty(cfg.resamplefs) && isempty(cfg.time),
  cfg.resamplefs = 256;
end

% select trials of interest
tmpcfg = keepfields(cfg, 'trials');
data   = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

if strcmp(cfg.sampleindex, 'yes') && isfield(data, 'sampleinfo')
  data.label{end+1} = 'sampleindex';
  for i=1:size(data.sampleinfo,1)
    % this works for one or more trials
    data.trial{i}(end+1,:) = data.sampleinfo(i,1):data.sampleinfo(i,2);
  end
elseif strcmp(cfg.sampleindex, 'yes')
  warning('no sampleinfo present, cannot add sampleindex as channel');
end

% sampleinfo, if present, becomes invalid because of the resampling
if isfield(data, 'sampleinfo'),
  data = rmfield(data, 'sampleinfo');
end

usefsample = ~isempty(cfg.resamplefs);
usetime    = ~isempty(cfg.time);

if usefsample && usetime
  error('you should either specify cfg.resamplefs or cfg.time')
end

% whether to use downsample() or resample()
usedownsample = 0;
if strcmp(cfg.resamplemethod, 'resample')
  usedownsample = 0;
elseif strcmp(cfg.resamplemethod, 'downsample')
  warning('using cfg.resamplemethod = ''downsample'', only use this if you have applied an anti-aliasing filter prior to downsampling!');
  usedownsample = 1;
else
  error('unknown resamplemethod ''%s''', cfg.resamplemethod);
end

% remember the original sampling frequency in the configuration
cfg.origfs = double(data.fsample);

if usefsample
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % resample based on new sampling frequency
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ntr = length(data.trial);

  ft_progress('init', cfg.feedback, 'resampling data');
  [fsorig, fsres] = rat(cfg.origfs./cfg.resamplefs);%account for non-integer fs
  cfg.resamplefs  = cfg.origfs.*(fsres./fsorig);%get new fs exact

  % make sure that the resampled time axes are aligned (this is to avoid
  % rounding errors in the time axes). this procedure relies on the
  % fact that resample assumes all data outside the data window to be zero
  % anyway. therefore, padding with zeros (to the left) before resampling
  % does not hurt
  firstsmp = zeros(ntr, 1);
  for itr = 1:ntr
    firstsmp(itr) = data.time{itr}(1);
  end
  minsmp = min(firstsmp);
  padsmp = round((firstsmp-minsmp).*cfg.origfs);

  nchan  = numel(data.label);
  if any(padsmp~=0)
    warning('not all of the trials have the same original time axis: to avoid rounding issues in the resampled time axes, data will be zero-padded to the left prior to resampling');
  end

  for itr = 1:ntr
    ft_progress(itr/ntr, 'resampling data in trial %d from %d\n', itr, ntr);
    if istrue(cfg.detrend)
      data.trial{itr} = ft_preproc_detrend(data.trial{itr});
    end

    % always remove the mean to avoid edge effects when there's a strong
    % offset, the cfg.demean option is dealt with below
    bsl             = nanmean(data.trial{itr},2);
    data.trial{itr} = data.trial{itr} - bsl(:,ones(1,size(data.trial{itr},2)));

    % pad the data with zeros to the left
    data.trial{itr} = [zeros(nchan, padsmp(itr))     data.trial{itr}];
    data.time{itr}  = [data.time{itr}(1)-(padsmp(itr):-1:1)./cfg.origfs data.time{itr}];

    % perform the resampling
    if usedownsample
      if mod(fsorig, fsres) ~= 0
        error('when using cfg.resamplemethod = ''downsample'', new sampling rate needs to be a proper divisor of original sampling rate');
      end

      if isa(data.trial{itr}, 'single')
        % temporary convert this trial to double precision
        data.trial{itr} = transpose(single(downsample(double(transpose(data.trial{itr})),fsorig/fsres)));
      else
        data.trial{itr} = transpose(downsample(transpose(data.trial{itr}),fsorig/fsres));
      end

    else % resample (default)
      if isa(data.trial{itr}, 'single')
        % temporary convert this trial to double precision
        data.trial{itr} = transpose(single(resample(double(transpose(data.trial{itr})),fsres,fsorig)));
      else
        data.trial{itr} = transpose(resample(transpose(data.trial{itr}),fsres,fsorig));
      end
    end

    % update the time axis
    nsmp = size(data.trial{itr},2);
    data.time{itr} = data.time{itr}(1) + (0:(nsmp-1))/cfg.resamplefs;

    %un-pad the data
    begindx         = ceil(cfg.resamplefs.*padsmp(itr)./cfg.origfs) + 1;
    data.time{itr}  = data.time{itr}(begindx:end);
    data.trial{itr} = data.trial{itr}(:, begindx:end);

    % add back the mean
    if ~strcmp(cfg.demean, 'yes')
      data.trial{itr} = data.trial{itr} + bsl(:,ones(1,numel(data.time{itr})));
    end

  end % for itr
  ft_progress('close');

  % specify the new sampling frequency in the output
  data.fsample = cfg.resamplefs;

elseif usetime
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % resample based on new time axes for each trial
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ntr = length(data.trial);

  ft_progress('init', cfg.feedback, 'resampling data');
  for itr = 1:ntr
    ft_progress(itr/ntr, 'resampling data in trial %d from %d\n', itr, ntr);

    if istrue(cfg.detrend)
      data.trial{itr} = ft_preproc_detrend(data.trial{itr});
    end

    % always remove the mean to avoid edge effects when there's a strong
    % offset, the cfg.demean option is dealt with below
    bsl             = nanmean(data.trial{itr},2);
    data.trial{itr} = data.trial{itr} - bsl(:,ones(1,size(data.trial{itr},2)));

    % perform the resampling
    if length(data.time{itr})>1,
      data.trial{itr} = interp1(data.time{itr}', data.trial{itr}', cfg.time{itr}', cfg.method)';
    else
      data.trial{itr} = repmat(data.trial{itr}, [1 length(cfg.time{itr}')]);
    end
    % update the time axis
    data.time{itr} = cfg.time{itr};

    % add back the mean
    if ~strcmp(cfg.demean, 'yes')
      data.trial{itr} = data.trial{itr} + bsl(:,ones(1,numel(data.time{itr})));
    end

  end % for itr
  ft_progress('close');

  % specify the new sampling frequency in the output
  t1 = cfg.time{1}(1);
  t2 = cfg.time{1}(2);
  data.fsample = 1/(t2-t1);

end % if usefsample or usetime

fprintf('original sampling rate = %d Hz\nnew sampling rate = %d Hz\n', cfg.origfs, data.fsample);

% convert back to input type if necessary
switch convert
  case 'timelock'
    data = ft_checkdata(data, 'datatype', 'timelock');
  otherwise
    % keep the output as it is
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data
