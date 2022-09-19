function [output] = volumethreshold(input, threshold, tissuelabel)

% VOLUMETHRESHOLD is a helper function for segmentations. It applies a
% relative threshold and subsequently looks for the largest connected part,
% thereby removing small blobs such as vitamine E capsules.
%
% See also VOLUMEFILLHOLES, VOLUMESMOOTH, VOLUMEPAD

if nargin<2 || isempty(threshold)
  threshold = 0;
end
if nargin<3 || isempty(tissuelabel)
  tissuelabel = 'volume';
end

% ensure that SPM is available, needed for spm_bwlabel
hasspm = ft_hastoolbox('spm8up', 3) || ft_hastoolbox('spm2', 1);

% mask by taking the negative of the segmentation, thus ensuring
% that no holes are within the compartment and do a two-pass
% approach to eliminate potential vitamin E capsules etc.

% ensure the input volume to be double precision, to allow for
% computations, unless it's already boolean
if ~isa(input, 'double') && ~isa(input, 'logical')
  input = double(input);
end

if ~isa(input, 'logical')
  if nargin==2, ft_error('if the input volume is not a boolean volume, you need to define a threshold value'); end
  if nargin==3, fprintf('thresholding %s at a relative threshold of %0.3f\n', tissuelabel, threshold); end
  output = double(input>(threshold*max(input(:))));
else
  % there is no reason to apply a threshold, but spm_bwlabel still needs a
  % double input for clustering
  output = double(input);
end

% cluster the connected tissue
[cluster, n] = spm_bwlabel(output, 6);

if n>1
  % it pays off to sort the cluster assignment if there are many clusters
  tmp = cluster(:);                       % convert to a vector
  tmp = tmp(tmp>0);                       % remove the zeros
  tmp = sort(tmp, 'ascend');              % sort according to cluster number
  m   = zeros(1,n);
  for k=1:n
    m(k) = sum(tmp==k);       % determine the voxel count for each cluster
    tmp  = tmp(m(k)+1:end);   % remove the last cluster that was counted
  end
  % select the tissue that has the most voxels belonging to it
  [m, i] = max(m);
  output = (cluster==i);
else
  % the output only contains a single cluster
  output = (cluster==1);
end
