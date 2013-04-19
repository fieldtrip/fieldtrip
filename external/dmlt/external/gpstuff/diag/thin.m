function x = thin(x,nburn,nthin,nlast)
%THIN  Delete burn-in and thin in MCMC-chains
%
%  Description
%    X = THIN(X,NBURN,NTHIN,NLAST) returns chain containing only
%    every NTHIN:th simulation sample starting from sample number
%    NBURN+1 and continuing to sample number NLAST.
%
%  See also
%    JOIN

% Copyright (c) 1999 Simo Särkkä
% Copyright (c) 2000,2010 Aki Vehtari
%
% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

if nargin < 4
  nlast = [];
end
if nargin < 3
  nthin = [];
end
if nargin < 2
  nburn = [];
end

[m,n]=size(x);
if isfield(x,'rstate')
  x=rmfield(x,'rstate');
end

if isstruct(x)
  if (m>1 | n>1)
    % array of structures
    for i=1:(m*n)
      x(i) = thin(x(i),nburn,nthin,nlast);
    end
  else
    % single structure
    names = fieldnames(x);
    for i=1:size(names,1)
      if isequal(names{i},'xtime')
        % Coxph model has ntime x 1 vector, which should be passed as is
        continue
      end
      value = getfield(x,names{i});
      if ~ischar(value) && (length(value) > 1 || isstruct(value))
	x = setfield(x,names{i},thin(value,nburn,nthin,nlast));
      elseif iscell(value)
	x = setfield(x,names{i},{thin(value{1},nburn,nthin,nlast)});
      end
    end
  end
elseif iscell(x)
  % cell array
  for i=1:(m*n)
    x{i} = thin(x{i},nburn,nthin,nlast);
  end
elseif m > 1
  % field array
  if isempty(nburn)
    nburn = 0;
  elseif (nburn < 0) | (nburn >= m)
    error('Illegal burn-in value');
  end
  if isempty(nthin)
    nthin = 1;
  elseif (nthin < 1) | (nthin > m)
    error('Illegal thinning value');
  end
  if isempty(nlast)
    nlast = m;
  elseif (nlast < 1) | (nlast > m)
    error('Illegal last index');
  end
  x = x((nburn+1):nthin:nlast,:);
end

end
