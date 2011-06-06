function [output, dimord] = ft_connectivity_bct(input, varargin)

% FT_CONNECTIVITY_BCT is a gateway function to the brain-connectivity-toolbox.

method = ft_getopt(varargin, 'method', '');
dimord = ft_getopt(varargin, 'dimord', '');

siz   = [size(input) 1];
input = reshape(input, [siz(1:2) prod(siz(3:end))]);

% check for binary or not
isbinary = true;
for k = 1:size(input,3)
  tmp = input(:,:,k);
  isbinary = all(ismember(tmp(:), [0 1]));
  if ~isbinary,
    break;
  end
end 

% check for directer or not
isdirected = true;
for k = 1:size(input,3)
  tmp = input(:,:,k);
  isdirected = all(all(tmp==tmp.'));
  if ~isdirected,
    break;
  end
end 

for k = 1:size(input, 3)
  switch method
  case 'assortativity'
  case 'betweenness'
  case 'breadthdist'
  case 'breadth'
  case 'charpath'
  case 'clustering_coef'

    % allocate memory
    if k==1 
      outsiz = size(input);
      outsiz(1) = []; 
      % remove one of the chan dimensions because the 
      % clustering coefficient is defined per channel
      % and not per channel pair
      output = zeros(outsiz);
      dimord = dimord(6:end);
    end

    if isbinary && isdirected
      output(:,k) = clustering_coef_bd(input(:,:,k));
    elseif isbinary && ~isdirected
      output(:,k) = clustering_coef_bu(input(:,:,k));
    elseif ~isbinary && isdirected
      output(:,k) = clustering_coef_wd(input(:,:,k));
    elseif ~isbinary && ~isdirected
      output(:,k) = clustering_coef_wu(input(:,:,k));
    end

  case 'degrees'
    
    % allocate memory
    if k==1
      outsiz = size(input);
      outsiz(1) = []; 
      % remove one of the chan dimensions because the 
      % degree is defined per channel
      % and not per channel pair
      output = zeros(outsiz);
      dimord = dimord(6:end);
    end

    if ~isbinary 
      warning_once('weights are not taken into account; graph is converted to binary values');
    end  
    
    if isdirected
      [in, out, output(:,k)] = degrees_dir(input(:,:,k));
      % fixme do something here
    elseif ~isdirected
      output(:,k) = degrees_und(input(:,:,k));
    end

  case 'density'
  case 'distance'
  case 'edge_betweenness'
  case 'efficiency'
  case 'modularity'
  case 'participation_coef'
  otherwise
    error('unsupported connectivity metric %s requested');
  end
end

%output = reshape(output, siz);
