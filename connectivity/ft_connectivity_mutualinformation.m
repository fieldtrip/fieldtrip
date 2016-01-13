function output = ft_connectivity_mutualinformation(input, varargin)

% computes mutual information using the information breakdown toolbox

ft_hastoolbox('ibtb', 1);

opts.nt     = [];
opts.method = 'dr';
opts.bias   = 'pt';

numbin  = ft_getopt(varargin, 'numbin', 10);
refindx = ft_getopt(varargin, 'refindx', 'all');
if ischar(refindx) && strcmp(refindx, 'all')
  refindx = (1:size(input,1))';
end

if length(refindx)~=numel(refindx)
  % could be channelcmb indexing
  error('channelcmb indexing is not, yet supported');
end

nsmp   = size(input,2);

output = zeros(numel(refindx), size(input,1))+nan;
for k = 1:numel(refindx)
  signal1 = input(refindx(k),:);
  
  % discretize signal1
  signal1 = binr(signal1, nsmp, numbin, 'eqpop');

  for m = setdiff(1:size(input,1),refindx(k))
    signal2 = input(m,:);
    
    % represent signal2 in bins according to signal1's discretization
    R = zeros(1,3,numbin);
    for j = 1:numbin
      nr         = signal1==j-1;
      opts.nt(j) = sum(nr);
      R(1, 1:opts.nt(j),j) = signal2(nr);
    end
    
    % discretize signal2 and compute mi
    R2 = binr(R, opts.nt', numbin, 'eqpop');
    output(k,m) = information(R2, opts, 'I'); % this computes mutual information
    
  end
end
