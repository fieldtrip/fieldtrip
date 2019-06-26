function [timelock] = ft_timelockbaseline(cfg, timelock)

% FT_TIMELOCKBASELINE performs baseline correction for ERF and ERP data
%
% Use as
%    [timelock] = ft_timelockbaseline(cfg, timelock)
% where the timelock data comes from FT_TIMELOCKANALYSIS and the
% configuration should contain
%   cfg.baseline     = [begin end] (default = 'no')
%   cfg.channel      = cell-array, see FT_CHANNELSELECTION
%   cfg.parameter    = field for which to apply baseline normalization, or
%                      cell-array of strings to specify multiple fields to normalize
%                      (default = 'avg')
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_TIMELOCKANALYSIS, FT_FREQBASELINE, FT_TIMELOCKGRANDAVERAGE

% Undocumented local options:
%   cfg.baselinewindow
%   cfg.previous
%   cfg.version

% Copyright (C) 2006, Robert Oostenveld
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
ft_preamble loadvar timelock
ft_preamble provenance timelock
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
timelock = ft_checkdata(timelock, 'datatype',...
  {'timelock+comp', 'timelock'}, 'feedback', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed', {'blc', 'demean'});
cfg = ft_checkconfig(cfg, 'renamed', {'blcwindow', 'baselinewindow'});
cfg = ft_checkconfig(cfg, 'forbidden', 'baselinetype');

% set the defaults
cfg.baseline  = ft_getopt(cfg, 'baseline',  'no');
cfg.channel   = ft_getopt(cfg, 'channel',   'all');
cfg.parameter = ft_getopt(cfg, 'parameter', '');

if isempty(cfg.parameter)
  if isfield(timelock, 'avg')
    cfg.parameter = 'avg';
  elseif strcmp(timelock.dimord, 'rpt_chan_time')
    cfg.parameter = 'trial';
  elseif strcmp(timelock.dimord, 'subj_chan_time')
    cfg.parameter = 'individual';
  end
end

% make sure cfg.parameter is a cell-array of strings
if (~isa(cfg.parameter, 'cell'))
  cfg.parameter = {cfg.parameter};
end

% the cfg.blc/blcwindow options are used in preprocessing and in
% ft_timelockanalysis (i.e. in private/preproc), hence make sure that
% they can also be used here for consistency
if isfield(cfg, 'baseline') && (isfield(cfg, 'demean') || isfield(cfg, 'baselinewindow'))
  ft_error('conflicting configuration options, you should use cfg.baseline');
elseif isfield(cfg, 'demean') && strcmp(cfg.demean, 'no')
  cfg.baseline = 'no';
  cfg = rmfield(cfg, 'demean');
  cfg = rmfield(cfg, 'baselinewindow');
elseif isfield(cfg, 'demean') && strcmp(cfg.demean, 'yes')
  cfg.baseline = cfg.baselinewindow;
  cfg = rmfield(cfg, 'demean');
  cfg = rmfield(cfg, 'baselinewindow');
end

if ischar(cfg.baseline)
  if strcmp(cfg.baseline, 'yes')
    % do correction on the whole time interval
    cfg.baseline = [-inf inf];
  elseif strcmp(cfg.baseline, 'all')
    % do correction on the whole time interval
    cfg.baseline = [-inf inf];
    % is input consistent?
  elseif strcmp(cfg.baseline, 'no')
    ft_warning('no baseline correction done');
    return;
  end
end

cfg.channel = ft_channelselection(cfg.channel, timelock.label);
chansel     = match_str(timelock.label, cfg.channel);

if ~(ischar(cfg.baseline) && strcmp(cfg.baseline, 'no'))
  % determine the time interval on which to apply baseline correction
  tbeg = nearest(timelock.time, cfg.baseline(1));
  tend = nearest(timelock.time, cfg.baseline(2));
  % update the configuration
  cfg.baseline(1) = timelock.time(tbeg);
  cfg.baseline(2) = timelock.time(tend);

   for k = 1:numel(cfg.parameter)
    par = cfg.parameter{k};

    % this if-statement is just there to give more specific text output
    if isequal(par, 'trial')
      fprintf('applying baseline correction on each individual trial\n');
      ntrial = size(timelock.(par),1);
      for i=1:ntrial
        timelock.(par)(i,chansel,:) = ft_preproc_baselinecorrect(shiftdim(timelock.(par)(i,chansel,:),1), tbeg, tend);
      end
    elseif isequal(par, 'individual')
      fprintf('applying baseline correction on each individual subject\n');
      nsubj = size(timelock.(par),1);
      for i=1:nsubj
        timelock.(par)(i,chansel,:) = ft_preproc_baselinecorrect(shiftdim(timelock.(par)(i,chansel,:),1), tbeg, tend);
      end
    else
      fprintf('applying baseline correction on %s\n', par);
      d = ndims(timelock.(par));
      if d == 3
        for i=1:size(timelock.(par),1)
          timelock.(par)(i,chansel,:) = ft_preproc_baselinecorrect(shiftdim(timelock.(par)(i,chansel,:),1), tbeg, tend);
        end
      elseif d == 2
        timelock.(par)(chansel,:) = ft_preproc_baselinecorrect(timelock.(par)(chansel,:), tbeg, tend);
      else
        ft_warning('Not doing anything - matrices up to only three dimensions are supported');
      end

    end

    if isfield(timelock, 'var')
      fprintf('baseline correction invalidates previous variance estimate, removing var\n');
      timelock = rmfield(timelock, 'var');
    end

    if isfield(timelock, 'cov')
      fprintf('baseline correction invalidates previous covariance estimate, removing cov\n');
      timelock = rmfield(timelock, 'cov');
    end
   end

end % ~strcmp(cfg.baseline, 'no')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output scaffolding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if numel(cfg.parameter)==1
  % convert from cell-array to string
  cfg.parameter = cfg.parameter{1};
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   timelock
ft_postamble provenance timelock
ft_postamble history    timelock
ft_postamble savevar    timelock
