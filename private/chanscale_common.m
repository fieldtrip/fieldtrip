function data = chanscale_common(cfg, data)

% CHANSCALE_COMMON applies a scaling to specific channel types
%
% Use as
%   data = chanscale_common(cfg, data)
% where the configuration contains
%   cfg.parameter
%
% For specific channel groups you can use
%   cfg.eegscale                = number, scaling to apply to the EEG channels prior to display
%   cfg.eogscale                = number, scaling to apply to the EOG channels prior to display
%   cfg.ecgscale                = number, scaling to apply to the ECG channels prior to display
%   cfg.emgscale                = number, scaling to apply to the EMG channels prior to display
%   cfg.megscale                = number, scaling to apply to the MEG channels prior to display
%   cfg.megrefscale             = number, scaling to apply to the MEG reference channels prior to display
%   cfg.magscale                = number, scaling to apply to the MEG magnetometer channels prior to display (in addition to the cfg.megscale factor)
%   cfg.gradscale               = number, scaling to apply to the MEG gradiometer channels prior to display (in addition to the cfg.megscale factor)
%   cfg.nirsscale               = number, scaling to apply to the NIRS channels prior to display
%
% For individual control off the scaling for all channels you can use
%   cfg.chanscale               = Nx1 vector with scaling factors, one per channel specified in cfg.channel
%
% For control over specific channels you can use
%   cfg.mychanscale             = number, scaling to apply to the channels specified in cfg.mychan
%   cfg.mychan                  = Nx1 cell-array with selection of channels

% Copyright (C) 2017-2019, Robert Oostenveld
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

cfg.eegscale    = ft_getopt(cfg, 'eegscale');
cfg.eogscale    = ft_getopt(cfg, 'eogscale');
cfg.ecgscale    = ft_getopt(cfg, 'ecgscale');
cfg.emgscale    = ft_getopt(cfg, 'emgscale');
cfg.megscale    = ft_getopt(cfg, 'megscale');
cfg.megrefscale = ft_getopt(cfg, 'megrefscale');
cfg.magscale    = ft_getopt(cfg, 'magscale');
cfg.gradscale   = ft_getopt(cfg, 'gradscale');
cfg.nirsscale   = ft_getopt(cfg, 'nirsscale');
cfg.chanscale   = ft_getopt(cfg, 'chanscale');
cfg.mychanscale = ft_getopt(cfg, 'mychanscale');
cfg.mychan      = ft_getopt(cfg, 'mychan');

% these should be column vectors
cfg.chanscale   = cfg.chanscale(:);
cfg.mychanscale = cfg.mychanscale(:);

if isfield(data, 'grad')
  % this helps to determine the different types of MEG channels
  senstype = ft_senstype(data.grad);
else
  senstype = [];
end

dimord = getdimord(data, cfg.parameter);

switch dimord
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case {'chan_time' 'chan_freq' 'chan_comp' 'chan_freq_time' 'chan_time_freq'}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get the data from the structure
    dat = data.(cfg.parameter);
    
    % apply scaling to selected channels, using wildcard to support subselection of channels
    if ~isempty(cfg.eegscale)
      chansel = match_str(data.label, ft_channelselection('EEG', data.label));
      ft_info('applying cfg.eegscale to %d channels', numel(chansel));
      dat(chansel,:,:) = dat(chansel,:,:) .* cfg.eegscale;
    end
    if ~isempty(cfg.eogscale)
      chansel = match_str(data.label, ft_channelselection('EOG', data.label));
      ft_info('applying cfg.eogscale to %d channels', numel(chansel));
      dat(chansel,:,:) = dat(chansel,:,:) .* cfg.eogscale;
    end
    if ~isempty(cfg.ecgscale)
      chansel = match_str(data.label, ft_channelselection('ECG', data.label));
      ft_info('applying cfg.ecgscale to %d channels', numel(chansel));
      dat(chansel,:,:) = dat(chansel,:,:) .* cfg.ecgscale;
    end
    if ~isempty(cfg.emgscale)
      chansel = match_str(data.label, ft_channelselection('EMG', data.label));
      ft_info('applying cfg.emgscale to %d channels', numel(chansel));
      dat(chansel,:,:) = dat(chansel,:,:) .* cfg.emgscale;
    end
    if ~isempty(cfg.megscale)
      chansel = match_str(data.label, ft_channelselection('MEG', data.label, senstype));
      ft_info('applying cfg.megscale to %d channels', numel(chansel));
      dat(chansel,:,:) = dat(chansel,:,:) .* cfg.megscale;
    end
    if ~isempty(cfg.megrefscale)
      chansel = match_str(data.label, ft_channelselection('MEGREF', data.label, senstype));
      ft_info('applying cfg.megrefscale to %d channels', numel(chansel));
      dat(chansel,:,:) = dat(chansel,:,:) .* cfg.megrefscale;
    end
    if ~isempty(cfg.magscale)
      chansel = match_str(data.label, ft_channelselection('MEGMAG', data.label, senstype));
      ft_info('applying cfg.magscale to %d channels', numel(chansel));
      dat(chansel,:,:) = dat(chansel,:,:) .* cfg.magscale;
    end
    if ~isempty(cfg.gradscale)
      chansel = match_str(data.label, ft_channelselection('MEGGRAD', data.label, senstype));
      ft_info('applying cfg.gradscale to %d channels', numel(chansel));
      dat(chansel,:,:) = dat(chansel,:,:) .* cfg.gradscale;
    end
    if ~isempty(cfg.nirsscale)
      chansel = match_str(data.label, ft_channelselection('NIRS', data.label));
      ft_info('applying cfg.nirsscale to %d channels', numel(chansel));
      dat(chansel,:,:) = dat(chansel,:,:) .* cfg.nirsscale;
    end
    if ~isempty(cfg.chanscale)
      chansel = match_str(data.label, ft_channelselection(cfg.channel, data.label));
      ft_info('applying cfg.chanscale to %d channels', numel(chansel));
      dat(chansel,:,:) = dat(chansel,:,:) .* repmat(cfg.chanscale,1,size(dat,2),size(dat,3));
    end
    if ~isempty(cfg.mychanscale)
      [chansel, scalesel] = match_str(data.label, cfg.mychan);
      ft_info('applying cfg.mychanscale to %d channels', numel(chansel));
      dat(chansel,:,:) = dat(chansel,:,:) .* repmat(cfg.mychanscale(scalesel),1,size(dat,2),size(dat,3));
    end
    
    % put the data back into the structure
    data.(cfg.parameter) = dat;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case '{rpt}_chan_time'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:numel(data.(cfg.parameter))
      % get the data from the structure
      dat = data.(cfg.parameter){i};
      
      % apply scaling to selected channels, using wildcard to support subselection of channels
      if ~isempty(cfg.eegscale)
        chansel = match_str(data.label, ft_channelselection('EEG', data.label));
        ft_info('applying cfg.eegscale to %d channels', numel(chansel));
        dat(chansel,:) = dat(chansel,:) .* cfg.eegscale;
      end
      if ~isempty(cfg.eogscale)
        chansel = match_str(data.label, ft_channelselection('EOG', data.label));
        ft_info('applying cfg.eogscale to %d channels', numel(chansel));
        dat(chansel,:) = dat(chansel,:) .* cfg.eogscale;
      end
      if ~isempty(cfg.ecgscale)
        chansel = match_str(data.label, ft_channelselection('ECG', data.label));
        ft_info('applying cfg.ecgscale to %d channels', numel(chansel));
        dat(chansel,:) = dat(chansel,:) .* cfg.ecgscale;
      end
      if ~isempty(cfg.emgscale)
        chansel = match_str(data.label, ft_channelselection('EMG', data.label));
        ft_info('applying cfg.emgscale to %d channels', numel(chansel));
        dat(chansel,:) = dat(chansel,:) .* cfg.emgscale;
      end
      if ~isempty(cfg.megscale)
        chansel = match_str(data.label, ft_channelselection('MEG', data.label, senstype));
        ft_info('applying cfg.megscale to %d channels', numel(chansel));
        dat(chansel,:) = dat(chansel,:) .* cfg.megscale;
      end
      if ~isempty(cfg.megrefscale)
        chansel = match_str(data.label, ft_channelselection('MEGREF', data.label, senstype));
        ft_info('applying cfg.megrefscale to %d channels', numel(chansel));
        dat(chansel,:) = dat(chansel,:) .* cfg.megrefscale;
      end
      if ~isempty(cfg.magscale)
        chansel = match_str(data.label, ft_channelselection('MEGMAG', data.label, senstype));
        ft_info('applying cfg.magscale to %d channels', numel(chansel));
        dat(chansel,:) = dat(chansel,:) .* cfg.magscale;
      end
      if ~isempty(cfg.gradscale)
        chansel = match_str(data.label, ft_channelselection('MEGGRAD', data.label, senstype));
        ft_info('applying cfg.gradscale to %d channels', numel(chansel));
        dat(chansel,:) = dat(chansel,:) .* cfg.gradscale;
      end
      if ~isempty(cfg.nirsscale)
        chansel = match_str(data.label, ft_channelselection('NIRS', data.label));
        ft_info('applying cfg.nirsscale to %d channels', numel(chansel));
        dat(chansel,:) = dat(chansel,:) .* cfg.nirsscale;
      end
      if ~isempty(cfg.chanscale)
        chansel = match_str(data.label, ft_channelselection(cfg.channel, data.label));
        ft_info('applying cfg.chanscale to %d channels', numel(chansel));
        dat(chansel,:) = dat(chansel,:) .* repmat(cfg.chanscale,1,size(dat,2));
      end
      if ~isempty(cfg.mychanscale)
        [chansel, scalesel] = match_str(data.label, cfg.mychan);
        ft_info('applying cfg.mychanscale to %d channels', numel(chansel));
        dat(chansel,:) = dat(chansel,:) .* repmat(cfg.mychanscale(scalesel),1,size(dat,2));
      end
      
      % put the data back into the structure
      data.(cfg.parameter){i} = dat;
      
      if i==1
        % only print information for the first trial
        ws = ft_info('off');
      end
      
    end % for each trial
    
    % revert to the normal information state
    ft_info(ws);
    
  otherwise
    if ~isempty(cfg.eegscale) || ...
        ~isempty(cfg.eogscale) || ...
        ~isempty(cfg.ecgscale) || ...
        ~isempty(cfg.emgscale) || ...
        ~isempty(cfg.megscale) || ...
        ~isempty(cfg.magscale) || ...
        ~isempty(cfg.gradscale) || ...
        ~isempty(cfg.nirsscale) || ...
        ~isempty(cfg.chanscale) || ...
        ~isempty(cfg.mychanscale)
      ft_error('unsupported dimord "%s"', dimord);
    end
end % switch dimord
