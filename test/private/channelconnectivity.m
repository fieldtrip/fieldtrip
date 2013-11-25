function [connectivity] = channelconnectivity(cfg, data)

% CHANNELCONNECTIVIY creates a NxN matrix that describes whether channels 
% are connected as neighbours
%
% See also FT_PREPARE_NEIGHBOURS

if (isfield(cfg, 'avgoverchan') && strcmp(cfg.avgoverchan, 'yes'))...
    || isempty(cfg.neighbours)
  connectivity = false(nchan,nchan);
else
  
  % check whether to use channels in the cfg or data
  if nargin < 2 && isfield(cfg, 'channel')
    chans = cfg.channel;
  elseif nargin == 2
    chans = data.label;
  else
    error('either the cfg needs to have both cfg.channel and cfg.neighbours, or a second (data) input argument needs to be specified');
  end
  
  nchan = length(chans);
  connectivity = false(nchan,nchan);
  
  % the original loop was the following, which is very slow:
%   for chan=1:length(cfg.neighbours)
%     [seld] = match_str(chans, cfg.neighbours(chan).label);
%     [seln] = match_str(chans, cfg.neighbours(chan).neighblabel);
%     if isempty(seld)
%       % this channel was not present in the data
%       continue;
%     else
%       % add the neighbours of this channel to the matrix
%       connectivity(seld, seln) = true;
%     end
%   end
  % the following code should be equivalent:
  
  
    
  % avoid using match_str inside the loop for performance reasons
  [sel1,sel2] = match_str(chans, {cfg.neighbours.label}');

  % make sure we only have the neighbours present in the data
  cfg.neighbours = cfg.neighbours(sel2);
  
  % make big list of all neighbour labels...
  allneighb = {cfg.neighbours.neighblabel};
  % ...store dimensionality...
  numAllNeighb = cellfun(@numel, allneighb);
  allneighb = cat(1, allneighb{:});
  % ...so that we can do with one call of match_str (with fullout=1)...
  selneighb = match_str(chans, allneighb, 1);
  % ...and then we can reshape the indices again to get the neighbours
  selneighb = mat2cell(selneighb, numAllNeighb);

  for k = 1:numel(sel1)
    % remove zeroes from selneighb (these correspond to channels absent)
    seln = selneighb{k};
    seln = seln(seln > 0);
    % add the neighbours of this channel to the matrix
    connectivity(sel1(k), seln) = true;
  end

end
