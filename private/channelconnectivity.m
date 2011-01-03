function [connectivity] = channelconnectivity(cfg)

% CHANNELCONNECTIVIY creates a NxN matrix that describes whether channels 
% are connected as neighbours
%
% See also FT_NEIGHBOURSELECTION

if strcmp(cfg.avgoverchan, 'yes')
  nchan = 1;
  connectivity = false(nchan,nchan);
else
  nchan = length(cfg.channel);
  connectivity = false(nchan,nchan);
  for chan=1:length(cfg.neighbours)
    [seld] = match_str(cfg.channel, cfg.neighbours{chan}.label);
    [seln] = match_str(cfg.channel, cfg.neighbours{chan}.neighblabel);
    if isempty(seld)
      % this channel was not present in the data
      continue;
    else
      % add the neighbours of this channel to the matrix
      connectivity(seld, seln) = true;
    end
  end
end
