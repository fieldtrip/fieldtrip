function [data] = ft_resampledata(cfg, data)

% FT_RESAMPLEDATA performs a resampling or downsampling of the data
%
% Use as
%   [data] = ft_resampledata(cfg, data)
%
% The data should be organised in a structure as obtained from
% the FT_PREPROCESSING function. The configuration should contain
%   cfg.resamplefs      = frequency at which the data will be resampled (default = 256 Hz)
%   cfg.detrend         = 'no' or 'yes', detrend the data prior to resampling (no default specified, see below)
%   cfg.demean          = 'no' or 'yes', whether to apply baseline correction (default = 'no')
%   cfg.baselinewindow  = [begin end] in seconds, the default is the complete trial (default = 'all')
%   cfg.feedback        = 'no', 'text', 'textbar', 'gui' (default = 'text')
%   cfg.trials          = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.sampleindex     = 'no' or 'yes', add a channel with the original sample indices (default = 'no')
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
% See also FT_PREPROCESSING, RESAMPLE, DOWNSAMPLE, INTERP1

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
cfg = ft_checkconfig(cfg, 'renamed', {'resamplemethod', 'method'});
cfg = ft_checkconfig(cfg, 'renamed', {'fsample', 'resamplefs'});

% set the defaults
cfg.resamplefs       = ft_getopt(cfg, 'resamplefs',      []);
cfg.time             = ft_getopt(cfg, 'time',            {});
cfg.factor           = ft_getopt(cfg, 'factor',          {});
cfg.detrend          = ft_getopt(cfg, 'detrend',         'no');
cfg.demean           = ft_getopt(cfg, 'demean',          'no');
cfg.baselinewindow   = ft_getopt(cfg, 'baselinewindow',  'all');
cfg.feedback         = ft_getopt(cfg, 'feedback',        'text');
cfg.trials           = ft_getopt(cfg, 'trials',          'all', 1);
cfg.method           = ft_getopt(cfg, 'method',          []);
cfg.sampleindex      = ft_getopt(cfg, 'sampleindex',     'no');

% store original datatype
convert = ft_datatype(data);

% check if the input data is valid for this function, this will convert it to raw if needed
data = ft_checkdata(data, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes');

if isempty(cfg.method) && ~isempty(cfg.time)
  % see INTERP1, shape-preserving piecewise cubic interpolation
  cfg.method = 'pchip';
elseif isempty(cfg.method)
  % see RESAMPLE
  cfg.method = 'resample';
end

usefsample = any(strcmp(cfg.method, {'resample', 'downsample', 'decimate', 'mean', 'median'}));
usetime    = ~usefsample;

% select trials of interest
tmpcfg = keepfields(cfg, {'trials', 'showcallinfo'});
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
  ft_warning('no sampleinfo present, cannot add sampleindex as channel');
end

% sampleinfo, if present, becomes invalid because of the resampling
if isfield(data, 'sampleinfo')
  data = rmfield(data, 'sampleinfo');
end

if usefsample && usetime
  ft_error('you should either specify cfg.resamplefs or cfg.time')
end

% remember the original sampling frequency in the configuration
cfg.origfs = double(data.fsample);

if usefsample
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % resample/downsample based on new sampling frequency
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
    ft_warning('not all of the trials have the same original time axis: to avoid rounding issues in the resampled time axes, data will be zero-padded to the left prior to resampling');
  end
  
  if any(strcmp(cfg.method, {'downsample', 'mean', 'median'}))
    ft_warning('using cfg.method = ''%s'', only use this if you have applied an anti-aliasing filter prior to downsampling!', cfg.method);
  end
  
  if any(strcmp(cfg.method, {'decimate', 'downsample', 'mean', 'median'}))
    if mod(fsorig, fsres) ~= 0
      ft_error('new sampling rate needs to be a proper divisor of original sampling rate');
    end
  end
  
  for itr = 1:ntr
    ft_progress(itr/ntr, 'resampling data in trial %d from %d\n', itr, ntr);
    if istrue(cfg.detrend)
      data.trial{itr} = ft_preproc_detrend(data.trial{itr});
    end
    
    % Compute mean based on a window (cfg.baselinewindow) or the whole
    % trial. This does not support nan values anymore. Fix in
    % ft_preproc_polyremoval.
    % previously --> bsl = nanmean(data.trial{itr},2);
    if ~strcmp(cfg.baselinewindow, 'all')
      begsample   = ceil(data.fsample*cfg.baselinewindow(1))+find(data.time{itr}==0);
      endsample   = floor(data.fsample*cfg.baselinewindow(2))-1+find(data.time{itr}==0);
      [~,bsl]     = ft_preproc_baselinecorrect(data.trial{itr}, begsample, endsample);
    else
      [~,bsl]     = ft_preproc_baselinecorrect(data.trial{itr});
    end
    % always remove the mean to avoid edge effects when there's a strong
    % offset, the cfg.demean option is dealt with below
    data.trial{itr} = data.trial{itr} - bsl(:,ones(1,size(data.trial{itr},2)));
    
    % pad the data with zeros to the left
    data.trial{itr} = [zeros(nchan, padsmp(itr))     data.trial{itr}];
    data.time{itr}  = [data.time{itr}(1)-(padsmp(itr):-1:1)./cfg.origfs data.time{itr}];
    
    % perform the resampling
    if strcmp(cfg.method, 'downsample')
      if isa(data.trial{itr}, 'single')
        % temporary convert this trial to double precision
        data.trial{itr} = transpose(single(downsample(double(transpose(data.trial{itr})),fsorig/fsres)));
      else
        data.trial{itr} = transpose(downsample(transpose(data.trial{itr}),fsorig/fsres));
      end
      
    elseif strcmp(cfg.method, 'resample')
      if isa(data.trial{itr}, 'single')
        % temporary convert this trial to double precision
        data.trial{itr} = transpose(single(resample(double(transpose(data.trial{itr})),fsres,fsorig)));
      else
        data.trial{itr} = transpose(resample(transpose(data.trial{itr}),fsres,fsorig));
      end
      
    elseif strcmp(cfg.method, 'decimate')
      if isa(data.trial{itr}, 'single')
        % temporary convert this trial to double precision
        data.trial{itr} = transpose(single(my_decimate(double(transpose(data.trial{itr})),fsorig/fsres)));
      else
        data.trial{itr} = transpose(my_decimate(transpose(data.trial{itr}),fsorig/fsres));
      end
      
    elseif strcmp(cfg.method, 'mean')
      if isa(data.trial{itr}, 'single')
        % temporary convert this trial to double precision
        data.trial{itr} = transpose(single(my_mean(double(transpose(data.trial{itr})),fsorig/fsres)));
      else
        data.trial{itr} = transpose(my_mean(transpose(data.trial{itr}),fsorig/fsres));
      end
      
    elseif strcmp(cfg.method, 'median')
      if isa(data.trial{itr}, 'single')
        % temporary convert this trial to double precision
        data.trial{itr} = transpose(single(my_median(double(transpose(data.trial{itr})),fsorig/fsres)));
      else
        data.trial{itr} = transpose(my_median(transpose(data.trial{itr}),fsorig/fsres));
      end
      
    else
      ft_error('unknown method ''%s''', cfg.method);
    end
    
    % update the time axis
    nsmp = size(data.trial{itr},2);
    origtime = data.time{itr};
    data.time{itr} = origtime(1) + (0:(nsmp-1))/cfg.resamplefs;
    % the first sample does not exactly remain at the same latency, but can be shifted by a sub-sample amount
    shift = mean(origtime) - mean(data.time{itr});
    data.time{itr} = data.time{itr} + shift;
    
    % un-pad the data
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
  % resample based on the specified new time axes for each trial
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ntr = length(data.trial);
  
  ft_progress('init', cfg.feedback, 'resampling data');
  for itr = 1:ntr
    ft_progress(itr/ntr, 'resampling data in trial %d from %d\n', itr, ntr);
    
    if istrue(cfg.detrend)
      data.trial{itr} = ft_preproc_detrend(data.trial{itr});
    end
    
    % Compute mean based on a window (cfg.baselinewindow) or the whole
    % trial. This does not support nan values anymore. Fix in
    % ft_preproc_polyremoval.
    % previously --> bsl = nanmean(data.trial{itr},2);
    if ~strcmp(cfg.baselinewindow, 'all')
      begsample   = ceil(data.fsample*cfg.baselinewindow(1))+find(data.time{itr}==0);
      endsample   = floor(data.fsample*cfg.baselinewindow(2))-1+find(data.time{itr}==0);
      [dum,bsl]   = ft_preproc_baselinecorrect(data.trial{itr}, begsample, endsample);
    else
      [dum,bsl]   = ft_preproc_baselinecorrect(data.trial{itr});
    end
    % always remove the mean to avoid edge effects when there's a strong
    % offset, the cfg.demean option is dealt with below
    data.trial{itr} = data.trial{itr} - bsl(:,ones(1,size(data.trial{itr},2)));
    
    % perform the resampling
    if length(data.time{itr})>1
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

ft_info('original sampling rate = %d Hz\nnew sampling rate = %d Hz\n', cfg.origfs, data.fsample);

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that decimates along the columns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = my_decimate(x, varargin)
[n, m] = size(x);
% decimate the first column
y = decimate(x(:,1), varargin{:});
if m>1
  % increase the size of the output matrix
  y(:,m) = 0;
  % decimate the subsequent columns
  for i=2:m
    y(:,i) = decimate(x(:,i), varargin{:});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that does a block-wise average along the columns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = my_mean(x, r)
[n, m] = size(x);
n = n - mod(n,r);
x = x(1:n,:);
x = reshape(x, [r n/r m]);
y = mean(x, 1);
y = reshape(y, [n/r m]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that does a block-wise median along the columns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = my_median(x, r)
[n, m] = size(x);
n = n - mod(n,r);
x = x(1:n,:);
x = reshape(x, [r n/r m]);
y = median(x, 1);
y = reshape(y, [n/r m]);
