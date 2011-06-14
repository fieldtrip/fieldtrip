function [pdc, pdcvar, n] = ft_connectivity_pdc(input, varargin)

hasjack  = keyval('hasjack',  varargin); if isempty(hasjack), hasjack = 0; end
powindx  = keyval('powindx',  varargin);
feedback = keyval('feedback', varargin); if isempty(feedback), feedback = 'none'; end
% FIXME build in proper documentation

% crossterms are described by chan_chan_therest
siz = [size(input) 1];
n   = siz(1);

outsum = zeros(siz(2:end));
outssq = zeros(siz(2:end));

% computing pdc is easiest on the inverse of the transfer function
pdim     = prod(siz(4:end));
tmpinput = reshape(input, [siz(1:3) pdim]);
ft_progress('init', feedback, 'inverting the transfer function...');
for k = 1:n
  ft_progress(k/n, 'inverting the transfer function for replicate %d from %d\n', k, n);
  tmp = reshape(tmpinput(k,:,:,:), [siz(2:3) pdim]);
  for m = 1:pdim
    tmp(:,:,m) = inv(tmp(:,:,m));
  end
  tmpinput(k,:,:,:) = tmp;
end
ft_progress('close');
input = reshape(tmpinput, siz);

ft_progress('init', feedback, 'computing metric...');
for j = 1:n
  ft_progress(j/n, 'computing metric for replicate %d from %d\n', j, n);
  invh   = reshape(input(j,:,:,:,:), siz(2:end));
  den    = sum(abs(invh).^2,1);
  tmppdc = abs(invh)./sqrt(repmat(den, [siz(2) 1 1 1 1]));
  %if ~isempty(cfg.submethod), tmppdc = baseline(tmppdc, cfg.submethod, baselineindx); end
  outsum = outsum + tmppdc;
  outssq = outssq + tmppdc.^2;
end
ft_progress('close');

pdc = outsum./n;

if n>1,
  if hasjack
    bias = (n-1).^2;
  else
    bias = 1;
  end
  pdcvar = bias*(outssq - (outsum.^2)./n)./(n - 1);
else
  pdcvar = [];
end

% this is to achieve the same convention for all connectivity metrics:
% row -> column
for k = 1:prod(siz(4:end))
  pdc(:,:,k)    = transpose(pdc(:,:,k));
  if ~isempty(pdcvar)
    pdcvar(:,:,k) = transpose(pdcvar(:,:,k));
  end
end
