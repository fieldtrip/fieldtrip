function [timelock] = ft_appendtimelock(cfg, varargin)

% FT_APPENDTIMELOCK concatenates multiple timelock (ERP/ERF) data
% structures that have been processed seperately. If the input data
% structures contain different channels, it will be concatenated along the
% channel direction. If the channels are identical in the input data
% structures, the data will be concatenated along the repetition dimension.
%
% Use as
%   combined = ft_appendtimelock(cfg, timelock1, timelock2, ...)
%
% The configuration can optionally contain
%   cfg.appenddim = String. The dimension to concatenate over (default:
%                   'auto').
%                   'chan' and 'rpt' possible (even if averaged data)
%   cfg.channel   = Nx1 cell-array with selection of channels (default = 'all'),
%                   see FT_CHANNELSELECTION for details
%
% See also FT_TIMELOCKANALYSIS, FT_APPENDDATA, FT_APPENDFREQ, FT_APPENDSOURCE

% Copyright (C) 2011, Robert Oostenveld
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
ft_preamble loadvar    varargin
ft_preamble provenance varargin
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

Ndata = length(varargin);

% check if the input data is valid for this function
for i=1:Ndata
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'timelock', 'feedback', 'yes', 'hassampleinfo', 'ifmakessense');
end

% set the defaults
cfg.channel    = ft_getopt(cfg, 'channel', 'all');
cfg.appenddim  = ft_getopt(cfg, 'appenddim', 'auto');
cfg.tolerance  = ft_getopt(cfg, 'tolerance',  1e-5);

% ensure that all inputs are sufficiently consistent
if ~checktime(varargin{:}, 'identical', cfg.tolerance);
  error('this function requires identical time axes for all input structures');
end

dimord = cell(1,Ndata);
for i=1:Ndata
  dimord{i} = varargin{i}.dimord;
end
dimordmatch = all(strcmp(dimord{1}, dimord));
if ~dimordmatch
  error('the dimords of the input data structures are not equal');
end

% determine over which dimension to append if not user-specified
switch cfg.appenddim
  case 'auto'
    [boolval, list] = checkchan(varargin{:}, 'unique');
    if boolval
      cfg.appenddim='chan';
    else
      cfg.appenddim='rpt';
    end
end

% start with the initial output structure
timelock        = [];
timelock.time   = varargin{1}.time;
ntime           = length(timelock.time);

switch cfg.appenddim
  case 'chan'
    for i=2:length(varargin) % check if any channels in common
      if ~isempty(ft_channelselection(varargin{1}.label,varargin{i}.label))
        error('Cannot concatenate over channels since some channels are in common');
      end
    end
    % no channels in common; append over channels
    nchan = zeros(size(varargin));
    for i=1:length(varargin)
      nchan(i) = length(varargin{i}.label);
    end
    chansel = cumsum([1 nchan]);
    timelock.label = cell(0,1);

    hascov        = isfield(varargin{1}, 'cov') && numel(size(varargin{1}.cov))==3;
    if hascov, warning('Concatenating over channels does not allow for keeping covariance information'); end

    if isfield(varargin{1}, 'trial')
      timelock.dimord = 'rpt_chan_time';

      % these don't make sense when concatenating the avg
      hastrialinfo  = isfield(varargin{1}, 'trialinfo');
      hassampleinfo = isfield(varargin{1}, 'sampleinfo');

      ntrial = size(varargin{1}.trial,1);


      timelock.trial = zeros(ntrial, sum(nchan), ntime);
      if hastrialinfo,  timelock.trialinfo = varargin{1}.trialinfo; end
      if hassampleinfo, timelock.sampleinfo = varargin{1}.sampleinfo; end

      for i=1:length(varargin)
        % copy the desired data into the output structure
        begchan = chansel(i);
        endchan = chansel(i+1)-1;
        timelock.trial(:,begchan:endchan,:) = varargin{i}.trial;
        timelock.label = [timelock.label; varargin{i}.label(:)];
      end % for varargin
      timelock.avg = permute(mean(timelock.trial,1),[2 3 1]);
      timelock.var = permute(var(timelock.trial,0,1),[2 3 1]);
      timelock.dimord = 'rpt_chan_time';

    else


      timelock.avg = zeros(sum(nchan), ntime);

      for i=1:length(varargin)
        % copy the desired data into the output structure
        begchan = chansel(i);
        endchan = chansel(i+1)-1;
        timelock.avg(begchan:endchan,:) = varargin{i}.avg;
        timelock.var(begchan:endchan,:) = varargin{i}.var;
        timelock.dof(begchan:endchan,:) = varargin{i}.dof;
        timelock.label = [timelock.label; varargin{i}.label(:)];
      end % for varargin

      % keep grad?  Should be ok for different channels of same dataset
      timelock.grad=varargin{1}.grad;
      timelock.dimord = 'chan_time';
    end


  case 'rpt' % append over trial or dataset dimension
    % select the channels that are in every dataset
    tmpcfg           = [];
    tmpcfg.channel   = cfg.channel;
    tmpcfg.tolerance = cfg.tolerance;
    [varargin{:}]    = ft_selectdata(tmpcfg, varargin{:});
    for i=1:Ndata
      [cfg_rolledback, varargin{i}] = rollback_provenance(cfg, varargin{i});
    end
    cfg = cfg_rolledback;

    timelock.label = varargin{1}.label;
    nchan          = numel(timelock.label);
    if nchan<1
      error('No channels in common');
    end

    hascov = isfield(varargin{1}, 'cov') && numel(size(varargin{1}.cov))==3;
    if isfield(varargin{1}, 'trial')
      % these don't make sense when concatenating the avg
      hastrialinfo  = isfield(varargin{1}, 'trialinfo');
      hassampleinfo = isfield(varargin{1}, 'sampleinfo');

      ntrial = zeros(size(varargin));
      for i=1:length(varargin)
        ntrial(i) = size(varargin{i}.trial, 1);
      end
      trialsel = cumsum([1 ntrial]);

      timelock.trial = zeros(sum(ntrial), nchan, ntime);
      if hastrialinfo,  timelock.trialinfo = zeros(sum(ntrial), size(varargin{1}.trialinfo,2)); end
      if hassampleinfo, timelock.sampleinfo = zeros(sum(ntrial), size(varargin{1}.sampleinfo,2)); end
      if hascov, timelock.cov = zeros(sum(ntrial), nchan, nchan); end

      for i=1:length(varargin)
        % copy the desired data into the output structure
        begtrial = trialsel(i);
        endtrial = trialsel(i+1)-1;
        chansel = match_str(cfg.channel, varargin{i}.label);
        timelock.trial(begtrial:endtrial,:,:) = varargin{i}.trial(:,chansel,:);
        if hastrialinfo,  timelock.trialinfo(begtrial:endtrial,:)   = varargin{i}.trialinfo(:,:); end
        if hassampleinfo, timelock.sampleinfo(begtrial:endtrial,:)  = varargin{i}.sampleinfo(:,:); end
        if hascov,        timelock.cov(begtrial:endtrial,:,:)       = varargin{i}.cov(:,chansel,chansel); end
      end % for varargin
      timelock.avg = permute(mean(timelock.trial,1),[2 3 1]);
      timelock.var = permute(var(timelock.trial,0,1),[2 3 1]);

    elseif isfield(varargin{1}, 'avg')

      ntrial = numel(varargin);
      timelock.trial = zeros(ntrial, nchan, ntime);
      if hascov, timelock.cov = zeros(sum(ntrial),nchan,nchan); end

      for i=1:length(varargin)
        % copy the desired data into the output structure
        chansel = match_str(cfg.channel, varargin{i}.label);
        timelock.trial(i,:,:) = varargin{i}.avg(chansel,:);
        if hascov, timelock.cov(i,:,:) = varargin{i}.cov(chansel,chansel); end
      end % for varargin
    end
    timelock.dimord = 'rpt_chan_time';
  otherwise
    error('it is not allowed to concatenate across dimension %s',cfg.appenddim);
end

% deal with the sensor information, if present
if isfield(varargin{1}, 'grad') || isfield(varargin{1}, 'elec')
  keepsensinfo = true;

  if isfield(varargin{1}, 'grad'), sensfield = 'grad'; end
  if isfield(varargin{1}, 'elec'), sensfield = 'elec'; end

  for k = 2:Ndata
    keepsensinfo = keepsensinfo && isequaln(varargin{1}.(sensfield), varargin{k}.(sensfield));
  end

  if keepsensinfo,
    timelock.(sensfield) = varargin{1}.(sensfield);
  end
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   varargin
ft_postamble provenance timelock
ft_postamble history    timelock
ft_postamble savevar    timelock
