function [c] = ft_connectivity_amplcorr_ortho(mom, varargin)

% Copyright (C) 2012 Jan-Mathijs Schoffelen

refindx = ft_getopt(varargin, 'refindx', 'all');
tapvec  = ft_getopt(varargin, 'tapvec',  ones(1,size(mom,2)));

if strcmp(refindx, 'all')
  refindx = 1:size(mom,1);
end

cmomnorm = conj(mom./abs(mom)); % only need to do conj() once

n        = size(mom,1);
ntap     = tapvec(1);
if ~all(tapvec==ntap)
  error('unequal number of tapers per observation is not yet supported');
end
if ntap>1
  error('more than one taper per observation is not yet supported');
end
tra  = zeros(size(mom,2), numel(tapvec));
for k = 1:numel(tapvec)
  tra((k-1)*ntap+(1:ntap), k) = 1./ntap;
end
powmom = (abs(mom).^2)*tra; % need only once
powmom = standardise(log10(powmom), 2);

c = zeros(n, numel(refindx)*2);
N = ones(n,1);
warning off;
for k = 1:numel(refindx)      
  indx     = refindx(k)
  ref      = mom(indx,:);
  crefnorm = conj(ref./abs(ref));

  % FIXME the following is probably not correct for ntap>1
  pow2 = (abs(imag(ref(N,:).*cmomnorm)).^2)*tra;
  pow2 = standardise(log10(pow2), 2);
  c1   = mean(powmom.*pow2, 2);
  pow1 = (abs(imag(mom.*crefnorm(N,:))).^2)*tra;
  pow2 = repmat((abs(ref).^2)*tra, [n 1]);
  pow1 = standardise(log10(pow1), 2);
  pow2 = standardise(log10(pow2), 2);
  c2   = mean(pow1.*pow2, 2);

  c(:,k) = c1;
  c(:,k+numel(refindx)) = c2;
end

