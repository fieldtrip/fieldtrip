function [data] = ft_resampledata(cfg, data);

% FT_RESAMPLEDATA performs a resampling or downsampling of the data
%
% Use as
%   [data] = ft_resampledata(cfg, data)
%
% The data should be organised in a structure as obtained from
% the FT_PREPROCESSING function. The configuration should contain
%   cfg.resamplefs = frequency at which the data will be resampled (default = 256 Hz)
%   cfg.detrend    = 'no' or 'yes', detrend the data prior to resampling (no default specified, see below)
%   cfg.blc        = 'no' or 'yes', baseline correct the data prior to resampling (default = 'no')
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
% See also FT_PREPROCESSING
%
% Undocumented local options:
%   cfg.inputfile  = one can specifiy preanalysed saved data as input
%   cfg.outputfile = one can specify output as file to save to disk

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

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

% set the defaults
if ~isfield(cfg, 'resamplefs'), cfg.resamplefs = [];   end
if ~isfield(cfg, 'time'),       cfg.time = {};         end
if ~isfield(cfg, 'detrend'),    cfg.detrend = [];      end  % no default to enforce people to consider backward compatibility problem, see below
if ~isfield(cfg, 'blc'),        cfg.blc = 'no';        end
if ~isfield(cfg, 'feedback'),   cfg.feedback = 'text'; end
if ~isfield(cfg, 'trials'),     cfg.trials = 'all';    end
if ~isfield(cfg, 'method'),     cfg.method = 'pchip';  end  % interpolation method
if ~isfield(cfg, 'inputfile'),  cfg.inputfile = [];    end
if ~isfield(cfg, 'outputfile'), cfg.outputfile = [];   end

% load optional given inputfile as data
hasdata = (nargin>1);
if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    data = loadvar(cfg.inputfile, 'data');
    hasdata = true;
  end
end

% check if the input data is valid for this function
data = checkdata(data, 'datatype', 'raw', 'feedback', 'yes');

if isempty(cfg.detrend)
  error('The previous default to apply detrending has been changed. Recommended is to apply a baseline correction instead of detrending. See the help of this function for more details.');
end

%set default resampling frequency
if isempty(cfg.resamplefs) && isempty(cfg.time),
  cfg.resamplefs = 256;
end

% this is needed if only a subset of trials is requested,
% and for an attempt to pass in the output a trl matrix
% which contains the indexing in the updated sampling rate
if isfield(data, 'cfg') % try to locate the trl in the nested configuration
  trl = findcfg(data.cfg, 'trl');
else
  trl = [];
end

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  fprintf('selecting %d trials\n', length(cfg.trials));
  data       = selectdata(data, 'rpt', cfg.trials);
  cfg.trlold = findcfg(data.cfg, 'trlold');
  cfg.trl    = findcfg(data.cfg, 'trl');
  trl        = cfg.trl;
end

usefsample = ~isempty(cfg.resamplefs);
usetime    = ~isempty(cfg.time);

if usefsample && usetime
  error('you should either specify cfg.resamplefs or cfg.time')
end

% remember the original sampling frequency in the configuration
cfg.origfs = data.fsample;

if usefsample
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % resample based on new sampling frequency
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ntr = length(data.trial);
  
  progress('init', cfg.feedback, 'resampling data');
  [fsorig, fsres] = rat(cfg.origfs./cfg.resamplefs);%account for non-integer fs
  cfg.resamplefs  = cfg.origfs.*(fsres./fsorig);%get new fs exact
  for itr = 1:ntr
    progress(itr/ntr, 'resampling data in trial %d from %d\n', itr, ntr);
    if strcmp(cfg.blc,'yes')
      data.trial{itr} = ft_preproc_baselinecorrect(data.trial{itr});
    end
    if strcmp(cfg.detrend,'yes')
      data.trial{itr} = ft_preproc_detrend(data.trial{itr});
    end
    % perform the resampling
    if isa(data.trial{itr}, 'single')
      % temporary convert this trial to double precision
      data.trial{itr} = single(resample(double(data.trial{itr})',fsres,fsorig))';
    else
      data.trial{itr} = resample(transpose(data.trial{itr}),fsres,fsorig)';
    end
    % update the time axis
    nsmp = size(data.trial{itr},2);
    data.time{itr} = data.time{itr}(1) + (0:(nsmp-1))/cfg.resamplefs;
  end % for itr
  progress('close');
  
  % specify the new sampling frequency in the output
  data.fsample = cfg.resamplefs;
  
elseif usetime
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % resample based on new time axes for each trial
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ntr = length(data.trial);
  
  progress('init', cfg.feedback, 'resampling data');
  for itr = 1:ntr
    progress(itr/ntr, 'resampling data in trial %d from %d\n', itr, ntr);
    if strcmp(cfg.blc,'yes')
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
  progress('close');
  
  % specify the new sampling frequency in the output
  t1 = cfg.time{1}(1);
  t2 = cfg.time{1}(2);
  data.fsample = 1/(t2-t1);
  
end % if usefsample or usetime

% update the trl matrix, which is needed for fetch_data and fetch_header (e.g. in ft_databrowser)
if ~isempty(trl)
  trlorig  = trl;
  offsindx = round((trl(:,1)-trl(:,3)).*(data.fsample./cfg.origfs));
  offs     = round(trl(:,3).*(data.fsample./cfg.origfs));
  nsmp     = cellfun('size',data.trial,2)';
  trl      = [offsindx+offs+1 offsindx+offs+nsmp offs];
  cfg.trl = trl;
  cfg.trlorig = trlorig;
end

fprintf('original sampling rate = %d Hz\nnew sampling rate = %d Hz\n', cfg.origfs, data.fsample);

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id$';

if hasdata && isfield(data, 'cfg')
  % remember the configuration details of the input data
  cfg.previous = data.cfg;
end
% remember the exact configuration details in the output
data.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', data); % use the variable name "data" in the output file
end

