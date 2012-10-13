function [data] = ft_resampledata(cfg, data)

% FT_RESAMPLEDATA performs a resampling or downsampling of the data
%
% Use as
%   [data] = ft_resampledata(cfg, data)
%
% The data should be organised in a structure as obtained from
% the FT_PREPROCESSING function. The configuration should contain
%   cfg.resamplefs = frequency at which the data will be resampled (default = 256 Hz)
%   cfg.detrend    = 'no' or 'yes', detrend the data prior to resampling (no default specified, see below)
%   cfg.demean     = 'no' or 'yes', baseline correct the data prior to resampling (default = 'no')
%   cfg.feedback   = 'no', 'text', 'textbar', 'gui' (default = 'text')
%   cfg.trials     = 'all' or a selection given as a 1xN vector (default = 'all')
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
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble provenance
ft_preamble trackconfig
ft_preamble loadvar data

% ft_checkdata is done further down

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed', {'blc', 'demean'});

% set the defaults
if ~isfield(cfg, 'resamplefs'), cfg.resamplefs = [];      end
if ~isfield(cfg, 'time'),       cfg.time       = {};      end
if ~isfield(cfg, 'detrend'),    cfg.detrend    = [];      end  % no default to enforce people to consider backward compatibility problem, see below
if ~isfield(cfg, 'demean'),     cfg.demean     = 'no';    end
if ~isfield(cfg, 'feedback'),   cfg.feedback   = 'text';  end
if ~isfield(cfg, 'trials'),     cfg.trials     = 'all';   end
if ~isfield(cfg, 'method'),     cfg.method     = 'pchip'; end  % interpolation method

% store original datatype
convert = ft_datatype(data);
  
% check if the input data is valid for this function, this will convert it to raw if needed
data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes');
  
if isempty(cfg.detrend)
  error('The previous default to apply detrending has been changed. Recommended is to apply a baseline correction instead of detrending. See the help of this function for more details.');
end

%set default resampling frequency
if isempty(cfg.resamplefs) && isempty(cfg.time),
  cfg.resamplefs = 256;
end

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  fprintf('selecting %d trials\n', length(cfg.trials));
  data       = ft_selectdata(data, 'rpt', cfg.trials);
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
    if strcmp(cfg.demean,'yes')
      data.trial{itr} = ft_preproc_baselinecorrect(data.trial{itr});
    end
    if strcmp(cfg.detrend,'yes')
      data.trial{itr} = ft_preproc_detrend(data.trial{itr});
    end

    % pad the data with zeros to the left
    data.trial{itr} = [zeros(nchan, padsmp(itr))     data.trial{itr}];
    data.time{itr}  = [data.time{itr}(1)-(padsmp(itr):-1:1)./cfg.origfs data.time{itr}];
    
    % perform the resampling
    if isa(data.trial{itr}, 'single')
      % temporary convert this trial to double precision
      data.trial{itr} = transpose(single(resample(double(transpose(data.trial{itr})),fsres,fsorig)));
    else
      data.trial{itr} = transpose(resample(transpose(data.trial{itr}),fsres,fsorig));
    end
    % update the time axis
    nsmp = size(data.trial{itr},2);
    data.time{itr} = data.time{itr}(1) + (0:(nsmp-1))/cfg.resamplefs;
  
    %un-pad the data
    begindx         = ceil(cfg.resamplefs.*padsmp(itr)./cfg.origfs) + 1;
    data.time{itr}  = data.time{itr}(begindx:end);
    data.trial{itr} = data.trial{itr}(:, begindx:end);
    
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
    if strcmp(cfg.demean,'yes')
      data.trial{itr} = ft_preproc_baselinecorrect(data.trial{itr});
    end
    if strcmp(cfg.detrend,'yes')
      data.trial{itr} = ft_preproc_detrend(data.trial{itr});
    end
    % perform the resampling
    if length(data.time{itr})>1,
      data.trial{itr} = interp1(data.time{itr}', data.trial{itr}', cfg.time{itr}', cfg.method)';
    else
      data.trial{itr} = repmat(data.trial{itr}, [1 length(cfg.time{itr}')]);
    end
    % update the time axis
    data.time{itr} = cfg.time{itr};
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
ft_postamble trackconfig
ft_postamble provenance
ft_postamble previous data
ft_postamble history data
ft_postamble savevar data
