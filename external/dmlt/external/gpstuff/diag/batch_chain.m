function x = batch(x,nburn,batchsize,nlast,fun)
%BATCH    Batch MCMC sample chain and evaluate mean/median of batches
%
%   X = BATCH(X, NBURN, BATCHSIZE, NLAST, FUN) Takes in sample chain X and 
%   omits NBURN samples from the beginning and divides it into batches of size 
%   BATCHSIZE until the NLAST'th sample. FUN is a function handle to handle the 
%   batches, for example @mean, which evaluates the mean of each batch after which 
%   a sample chain X is returned containing the mean values of each batch. User can
%   give any function for batch as a function handle.
%
%   The default values are BATCHSIZE = 2, NLAST = SIZE(X,1) and FUN = @mean.
%
%   See also
%     THIN

% Copyright (c) 1999 Simo Särkkä
% Copyright (c) 2000 Aki Vehtari
% Copyright (c) 2004 Jarno Vanhatalo
%
% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

if nargin < 5
  fun = @mean;
end

if nargin < 4
  nlast = [];
end

if nargin < 3
  batchsize = [];
end

[m,n]=size(x);
if isfield(x,'rstate')
  x=rmfield(x,'rstate');
end

if isstruct(x)
  if (m>1 | n>1)
    % array of structures
    for i=1:(m*n)
      x(i) = batch(x(i),nburn,batchsize,nlast,fun);
    end
  else
    % single structure
    names = fieldnames(x);
    for i=1:size(names,1)
      value = getfield(x,names{i});
      if length(value) > 1
	x = setfield(x,names{i},batch(value,nburn,batchsize,nlast,fun));
      elseif iscell(value)
	x = setfield(x,names{i},{batch(value{1},nburn,batchsize,nlast,fun)});
      end
    end
  end
elseif iscell(x)
  % cell array
  for i=1:(m*n)
    x{i} = batch(x{i},nburn,batchsize,nlast,fun);
  end
elseif m > 1
  % field array
  if isempty(nlast)
    nlast = m;
  elseif (nlast < 0) | (nlast >= m)
    error('Illegal nlast value');
  end
  if isempty(nburn)
    nburn = 0;
  elseif (nburn < 0) | (nburn >= nlast)
    error('Illegal burn-in value');
  end
  if isempty(batchsize)
    batchsize = 1;
  elseif (batchsize <= 1) | (batchsize > m)
    error('Illegal batchsize value');
  end
  if isempty(fun)
    fun = 1;
  end
  
  B = [];
  for i=1:ceil((nlast-nburn)/batchsize)
    if (nburn+i*batchsize) < nlast
      B(i,:) = feval(fun,(x((nburn+1+(i-1)*batchsize):(nburn+i*batchsize),:)));
    else
      B(i,:) = feval(fun,(x((nburn+1+(i-1)*batchsize):nlast,:)));
    end
  end
  x=B;  
end
