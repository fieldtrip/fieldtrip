function [connectivity] = channelconnectivity(cfg, data)

% CHANNELCONNECTIVIY creates a NxN matrix that describes whether channels
% are connected as neighbours
%
% See also FT_PREPARE_NEIGHBOURS

if (isfield(cfg, 'avgoverchan') && strcmp(cfg.avgoverchan, 'yes')) || isempty(cfg.neighbours)
  if nargin < 2
    nchan = numel(cfg.channel);
  else
    nchan = numel(data.label);
  end
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
  
  nchan        = length(chans);
  connectivity = false(nchan,nchan);
  neighbours   = struct(cfg.neighbours);
  
  try
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
    % the above is still used in case the fast version fails (catch below)
    % the following code should be equivalent:
    
    if numel(neighbours) < nchan/4
      error('channelconnectivity only works when at least 1/4th of all channels has neighbours defined (or with neighbours = [])');
      % FIXME the above error is true (or maybe with <1/5th) but I don't
      % understand why, ES 25-nov-2013
      % the above error will cause the code to fall through to the catch
      % below, so no big worries
    end
    
    % avoid using match_str inside the loop for performance reasons
    [sel1,sel2] = match_str(chans, {neighbours.label}');
    
    % make sure we only have the neighbours present in the data
    neighbours = neighbours(sel2);
    
    % make big list of all neighbour labels...
    allneighb = {neighbours.neighblabel};
    
    % determine whether neighbours are stored as row or column vector
    siz1 = cellfun(@size, allneighb, repmat({1},size(allneighb)));
    siz2 = cellfun(@size, allneighb, repmat({2},size(allneighb)));
    if sum(siz2) > sum(siz1)
      catdim = 2;
    else
      catdim = 1;
    end
    
    % ...store dimensionality...
    numAllNeighb = cellfun(@numel, allneighb);
    allneighb = cat(catdim, allneighb{:});
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
    
  catch
    % the fast version failed, use slow instead
    
    for chan=1:length(neighbours)
      [seld] = match_str(chans, neighbours(chan).label);
      [seln] = match_str(chans, neighbours(chan).neighblabel);
      if isempty(seld)
        % this channel was not present in the data
        continue;
      else
        % add the neighbours of this channel to the matrix
        connectivity(seld, seln) = true;
      end
    end
    
  end
end
