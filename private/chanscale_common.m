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
%   cfg.gradscale               = number, scaling to apply to the MEG gradiometer channels prior to display (in addition to the cfg.megscale factor)
%   cfg.magscale                = number, scaling to apply to the MEG magnetometer channels prior to display (in addition to the cfg.megscale factor)
%
% For individual control off the scaling for all channels you can use
%   cfg.chanscale               = Nx1 vector with scaling factors, one per channel specified in cfg.channel
%
% For control over specific channels you can use
%   cfg.mychanscale             = number, scaling to apply to the channels specified in cfg.mychan
%   cfg.mychan                  = Nx1 cell-array with selection of channels

cfg.eegscale    = ft_getopt(cfg, 'eegscale');
cfg.eogscale    = ft_getopt(cfg, 'eogscale');
cfg.ecgscale    = ft_getopt(cfg, 'ecgscale');
cfg.emgscale    = ft_getopt(cfg, 'emgscale');
cfg.megscale    = ft_getopt(cfg, 'megscale');
cfg.gradscale   = ft_getopt(cfg, 'gradscale');
cfg.magscale    = ft_getopt(cfg, 'magscale');
cfg.chanscale   = ft_getopt(cfg, 'chanscale');
cfg.mychanscale = ft_getopt(cfg, 'mychanscale');
cfg.mychan      = ft_getopt(cfg, 'mychan');

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
      dat(chansel,:,:) = dat(chansel,:,:) .* cfg.eegscale;
    end
    if ~isempty(cfg.eogscale)
      chansel = match_str(data.label, ft_channelselection('EOG', data.label));
      dat(chansel,:,:) = dat(chansel,:,:) .* cfg.eogscale;
    end
    if ~isempty(cfg.ecgscale)
      chansel = match_str(data.label, ft_channelselection('ECG', data.label));
      dat(chansel,:,:) = dat(chansel,:,:) .* cfg.ecgscale;
    end
    if ~isempty(cfg.emgscale)
      chansel = match_str(data.label, ft_channelselection('EMG', data.label));
      dat(chansel,:,:) = dat(chansel,:,:) .* cfg.emgscale;
    end
    if ~isempty(cfg.megscale)
      type = opt.hdr.grad.type;
      chansel = match_str(data.label, ft_channelselection('MEG', data.label, type));
      dat(chansel,:,:) = dat(chansel,:,:) .* cfg.megscale;
    end
    if ~isempty(cfg.magscale)
      chansel = match_str(data.label, ft_channelselection('MEGMAG', data.label));
      dat(chansel,:,:) = dat(chansel,:,:) .* cfg.magscale;
    end
    if ~isempty(cfg.gradscale)
      chansel = match_str(data.label, ft_channelselection('MEGGRAD', data.label));
      dat(chansel,:,:) = dat(chansel,:,:) .* cfg.gradscale;
    end
    if ~isempty(cfg.chanscale)
      chansel = match_str(data.label, ft_channelselection(cfg.channel, data.label));
      dat(chansel,:,:) = dat(chansel,:,:) .* repmat(cfg.chanscale(:),1,size(dat,2),size(dat,3));
    end
    if ~isempty(cfg.mychanscale)
      chansel = match_str(data.label, ft_channelselection(cfg.mychan, data.label));
      dat(chansel,:,:) = dat(chansel,:,:) .* repmat(cfg.mychanscale(:),1,size(dat,2),size(dat,3));
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
        dat(chansel,:) = dat(chansel,:) .* cfg.eegscale;
      end
      if ~isempty(cfg.eogscale)
        chansel = match_str(data.label, ft_channelselection('EOG', data.label));
        dat(chansel,:) = dat(chansel,:) .* cfg.eogscale;
      end
      if ~isempty(cfg.ecgscale)
        chansel = match_str(data.label, ft_channelselection('ECG', data.label));
        dat(chansel,:) = dat(chansel,:) .* cfg.ecgscale;
      end
      if ~isempty(cfg.emgscale)
        chansel = match_str(data.label, ft_channelselection('EMG', data.label));
        dat(chansel,:) = dat(chansel,:) .* cfg.emgscale;
      end
      if ~isempty(cfg.megscale)
        type = opt.hdr.grad.type;
        chansel = match_str(data.label, ft_channelselection('MEG', data.label, type));
        dat(chansel,:) = dat(chansel,:) .* cfg.megscale;
      end
      if ~isempty(cfg.magscale)
        chansel = match_str(data.label, ft_channelselection('MEGMAG', data.label));
        dat(chansel,:) = dat(chansel,:) .* cfg.magscale;
      end
      if ~isempty(cfg.gradscale)
        chansel = match_str(data.label, ft_channelselection('MEGGRAD', data.label));
        dat(chansel,:) = dat(chansel,:) .* cfg.gradscale;
      end
      if ~isempty(cfg.chanscale)
        chansel = match_str(data.label, ft_channelselection(cfg.channel, data.label));
        dat(chansel,:) = dat(chansel,:) .* repmat(cfg.chanscale(:),1,size(dat,2));
      end
      if ~isempty(cfg.mychanscale)
        chansel = match_str(data.label, ft_channelselection(cfg.mychan, data.label));
        dat(chansel,:) = dat(chansel,:) .* repmat(cfg.mychanscale(:),1,size(dat,2));
      end
      
      % put the data back into the structure
      data.(cfg.parameter){i} = dat;
    end % for each trial
    
  otherwise
    if ~isempty(cfg.eegscale) || ...
        ~isempty(cfg.eogscale) || ...
        ~isempty(cfg.ecgscale) || ...
        ~isempty(cfg.emgscale) || ...
        ~isempty(cfg.megscale) || ...
        ~isempty(cfg.magscale) || ...
        ~isempty(cfg.gradscale) || ...
        ~isempty(cfg.chanscale) || ...
        ~isempty(cfg.mychanscale)
      ft_error('unsupported dimord "%s"', dimord);
    end
end % switch dimord
