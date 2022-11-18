function [data] = ft_resampledata(cfg, data)

% FT_RESAMPLEDATA performs a resampling or downsampling of the data
%
% Use as
%   [data] = ft_resampledata(cfg, data)
%
% The data should be organised in a structure as obtained from the FT_PREPROCESSING
% function. The configuration should contain
%   cfg.resamplefs      = frequency at which the data will be resampled
%   cfg.detrend         = 'no' or 'yes', detrend the data prior to resampling (no default specified, see below)
%   cfg.demean          = 'no' or 'yes', whether to apply baseline correction (default = 'no')
%   cfg.baselinewindow  = [begin end] in seconds, the default is the complete trial (default = 'all')
%   cfg.feedback        = 'no', 'text', 'textbar', 'gui' (default = 'text')
%   cfg.trials          = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.sampleindex     = 'no' or 'yes', add a channel with the original sample indices (default = 'no')
%
% Instead of specifying cfg.resamplefs, you can also specify a time axis on which you
% want the data to be resampled. This is useful for merging data from two acquisition
% devices, after resampledata you can call FT_APPENDDATA to concatenate the channels
% from the different acquisition devices.
%   cfg.time        = cell-array with one time axis per trial (i.e. from another dataset)
%   cfg.method      = interpolation method, see INTERP1 (default = 'pchip')
%   cfg.extrapval   = extrapolation behaviour, scalar value or 'extrap' (default = as in INTERP1)
%
% When you specify cfg.method='resample' an implicit anti-aliasing low pass filter is applied prior
% to the resampling. You can also explicitly specify an anti-aliasing low pass filter. This is adviced
% when downsampling using any other method than 'resample', but also when strong noise components are
% present just above the new Nyquist frequency.
%   cfg.lpfilter    = 'yes' or 'no' (default = 'no')
%   cfg.lpfreq      = scalar value for low pass frequency (there is no default, so needs to be always specified)
%   cfg.lpfilttype  = string, filter type (default is set in ft_preproc_lowpassfilter)
%   cfg.lpfiltord   = scalar, filter order (default is set in ft_preproc_lowpassfilter)
%
% More documentation about anti-alias filtering can be found in this <a href="matlab:
% web('https://www.fieldtriptoolbox.org/faq/resampling_lowpassfilter')">FAQ</a> on the FieldTrip website.
%
% Previously this function used to detrend the data by default to avoid edge artifacts fue to the
% anti-aliassing filter. Detrending is fine for removing slow drifts in data prior to frequency analysis,
% but not recommended if you want to look at ERPs or ERFs. Therefore the old default value 'yes' has been
% removed; you now explicitly have to specify whether you want to detrend or not.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_PREPROCESSING, FT_APPENDDATA, FT_PREPROC_LOWPASSFILTER, ESAMPLE, DOWNSAMPLE, INTERP1

% Copyright (C) 2003-2006, FC Donders Centre, Markus Siegel
% Copyright (C) 2004-2022, FC Donders Centre, Robert Oostenveld
% Copyright (C) 2022, DCCN, Jan-Mathijs Schoffelen
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

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% ft_checkdata is done further down

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden', {'trial'}); % prevent accidental typos, see issue 1729
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
cfg.extrapval        = ft_getopt(cfg, 'extrapval',       []);
cfg.lpfilter         = ft_getopt(cfg, 'lpfilter');

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
tmpcfg = keepfields(cfg, {'trials', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
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

% set this to nan, it will be updated later on
data.fsample = nan;

if isempty(cfg.lpfilter), cfg.lpfilter = 'no'; end
dolpfilt = istrue(cfg.lpfilter);
if dolpfilt
  cfg.lpfilttype = ft_getopt(cfg, 'lpfilttype');
  cfg.lpfiltord  = ft_getopt(cfg, 'lpfiltord');
  cfg            = ft_checkconfig(cfg, 'required', 'lpfreq');
end

if usefsample
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % resample/downsample based on new sampling frequency
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ntr = length(data.trial);
  nchan  = numel(data.label);

  ft_progress('init', cfg.feedback, 'resampling data');
  [fsorig, fsres] = rat(cfg.origfs./cfg.resamplefs); %account for non-integer fs
  cfg.resamplefs  = cfg.origfs.*(fsres./fsorig); %get new fs exact

  % make sure that the resampled time axes are aligned (this is to avoid rounding
  % errors in the time axes). this procedure relies on the fact that resample assumes
  % all data outside the data window to be zero anyway. therefore, padding with zeros
  % (to the left and right) before resampling does not hurt
  begsample = zeros(ntr, 1);
  endsample = zeros(ntr, 1);
  for itr = 1:ntr
    begsample(itr) = round(cfg.origfs * data.time{itr}(1));
    endsample(itr) = round(cfg.origfs * data.time{itr}(end));
  end
  begpad = begsample-min(begsample);
  endpad = max(endsample)-endsample;

  if any(begpad~=0) || any(endpad~=0)
    ft_warning('not all trials have the same time axis; data will be zero-padded prior to resampling to avoid rounding issues in the resampled time axes');
  end

  if any(strcmp(cfg.method, {'downsample', 'mean', 'median'}))
    ft_warning('using cfg.method = ''%s''; only use this if you have applied an anti-aliasing filter prior to downsampling!', cfg.method);
  end

  if any(strcmp(cfg.method, {'decimate', 'downsample', 'mean', 'median'}))
    if mod(fsorig, fsres) ~= 0
      ft_error('the new sampling rate needs to be an integer division of the original sampling rate');
    end
  end

  for itr = 1:ntr
    ft_progress(itr/ntr, 'resampling data in trial %d from %d\n', itr, ntr);

    olddat = data.trial{itr};
    oldtim = data.time{itr};

    % detrending is in general not recommended
    if istrue(cfg.detrend)
      if ~strcmp(cfg.baselinewindow, 'all')
        olddat = ft_preproc_detrend(olddat, nearest(oldtim, cfg.baselinewindow(1)), nearest(oldtim, cfg.baselinewindow(2)));
      else
        olddat = ft_preproc_detrend(olddat);
      end
    end

    % remove the mean to avoid edge effects when there's a strong offset, the cfg.demean option is dealt with below
    if ~strcmp(cfg.baselinewindow, 'all')
      [olddat, bsl] = ft_preproc_baselinecorrect(olddat, nearest(oldtim, cfg.baselinewindow(1)), nearest(oldtim, cfg.baselinewindow(2)));
    else
      [olddat, bsl] = ft_preproc_baselinecorrect(olddat);
    end

    if istrue(cfg.lpfilter)
      olddat = ft_preproc_lowpassfilter(olddat, cfg.origfs, cfg.lpfreq, cfg.lpfiltord, cfg.lpfilttype);
    end

    % pad the data with zeros on both sides
    olddat = [zeros(nchan, begpad(itr)) olddat zeros(nchan, endpad(itr))];
    oldtim = ((begsample(itr)-begpad(itr)):(endsample(itr)+endpad(itr))) / cfg.origfs;

    % perform the resampling
    if strcmp(cfg.method, 'downsample')
      if isa(olddat, 'single')
        % temporary convert this trial to double precision
        newdat = transpose(single(downsample(double(transpose(olddat)),fsorig/fsres)));
      else
        newdat = transpose(downsample(transpose(olddat),fsorig/fsres));
      end

    elseif strcmp(cfg.method, 'resample')
      if isa(olddat, 'single')
        % temporary convert this trial to double precision
        newdat = transpose(single(resample(double(transpose(olddat)),fsres,fsorig)));
      else
        newdat = transpose(resample(transpose(olddat),fsres,fsorig));
      end

    elseif strcmp(cfg.method, 'decimate')
      if isa(olddat, 'single')
        % temporary convert this trial to double precision
        newdat = transpose(single(my_decimate(double(transpose(olddat)),fsorig/fsres)));
      else
        newdat = transpose(my_decimate(transpose(olddat),fsorig/fsres));
      end

    elseif strcmp(cfg.method, 'mean')
      if isa(olddat, 'single')
        % temporary convert this trial to double precision
        newdat = transpose(single(my_mean(double(transpose(olddat)),fsorig/fsres)));
      else
        newdat = transpose(my_mean(transpose(olddat),fsorig/fsres));
      end

    elseif strcmp(cfg.method, 'median')
      if isa(olddat, 'single')
        % temporary convert this trial to double precision
        newdat = transpose(single(my_median(double(transpose(olddat)), fsorig/fsres)));
      else
        newdat = transpose(my_median(transpose(olddat), fsorig/fsres));
      end

    else
      ft_error('unknown method ''%s''', cfg.method);
    end

    % add back the mean
    if ~strcmp(cfg.demean, 'yes')
      nsmp   = size(newdat,2);
      newdat = newdat + bsl(:,ones(1,nsmp));
    end

    % compute the new time axis, assuming that it starts at the same time
    nsmp   = size(newdat,2);
    newtim = (0:(nsmp-1))/cfg.resamplefs;

    % the middle of the time bin represented by the first samples are not aligned
    % the new time axis can be shifted by a sub-sample amount
    shift  = mean(oldtim) - mean(newtim);
    newtim = newtim + shift;

    if begpad(itr)>0 || endpad(itr)>0
      % un-pad the data
      sel = (1+round(begpad(itr)*cfg.resamplefs/cfg.origfs)):(length(newtim)-round(endpad(itr)*cfg.resamplefs/cfg.origfs));
      newtim = newtim(   sel);
      newdat = newdat(:, sel);
    end

    data.time{itr}  = newtim;
    data.trial{itr} = newdat;

  end % for itr
  ft_progress('close');

  % specify the new sampling frequency in the output
  data.fsample = cfg.resamplefs;

elseif usetime
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % resample based on the specified new time axes for each trial
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if isempty(cfg.extrapval)
    if strcmp(cfg.method, 'spline') || strcmp(cfg.method, 'pchip')
      cfg.extrapval = 'extrap';
    else
      cfg.extrapval = nan;
    end
  end

  ntr = length(data.trial);

  ft_progress('init', cfg.feedback, 'resampling data');
  for itr = 1:ntr
    ft_progress(itr/ntr, 'resampling data in trial %d from %d\n', itr, ntr);

    olddat = data.trial{itr};
    oldtim = data.time{itr};

    % detrending is in general not recommended
    if istrue(cfg.detrend)
      if ~strcmp(cfg.baselinewindow, 'all')
        olddat = ft_preproc_detrend(olddat, nearest(oldtim, cfg.baselinewindow(1)), nearest(oldtim, cfg.baselinewindow(2)));
      else
        olddat = ft_preproc_detrend(olddat);
      end
    end

    % always remove the mean to avoid edge effects when there's a strong offset, the cfg.demean option is dealt with below
    if ~strcmp(cfg.baselinewindow, 'all')
      [olddat, bsl] = ft_preproc_baselinecorrect(olddat, nearest(oldtim, cfg.baselinewindow(1)), nearest(oldtim, cfg.baselinewindow(2)));
    else
      [olddat, bsl] = ft_preproc_baselinecorrect(olddat);
    end

    if istrue(cfg.lpfilter)
      olddat = ft_preproc_lowpassfilter(olddat, cfg.origfs, cfg.lpfreq, cfg.lpfiltord, cfg.lpfilttype);
    end

    % perform the resampling
    newtim = cfg.time{itr};
    if length(oldtim)>1
      newdat = interp1(oldtim', olddat', newtim', cfg.method, cfg.extrapval)';
    else
      newdat = repmat(olddat, [1 numel(newtim)]);
    end

    % add back the mean
    if ~strcmp(cfg.demean, 'yes')
      nsmp   = size(newdat, 2);
      newdat = newdat + bsl(:,ones(1,nsmp));
    end

    data.trial{itr} = newdat;
    data.time{itr}  = newtim;

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
