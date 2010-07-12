function [c, v, n] = ft_connectivity_psi(cfg, input, hasrpt, hasjack)

if nargin==2,
  hasrpt   = 0;
  hasjack  = 0;
elseif nargin==3,
  hasjack  = 0;
end

if (length(strfind(cfg.dimord, 'chan'))~=2 || length(strfind(cfg.dimord, 'pos'))>0) && isfield(cfg, 'powindx') && ~isempty(cfg.powindx),
  %crossterms are not described with chan_chan_therest, but are linearly indexed
  
  siz = size(input);
  if ~hasrpt,
    siz   = [1 siz];
    input = reshape(input, siz);
  end
  
  outsum = zeros(siz(2:end));
  outssq = zeros(siz(2:end));
  pvec   = [2 setdiff(1:numel(siz),2)];
  
  progress('init', cfg.feedback, 'computing metric...');
  %first compute coherency and then phaseslopeindex
  for j = 1:siz(1)
    progress(j/siz(1), 'computing metric for replicate %d from %d\n', j, siz(1));
    c      = reshape(input(j,:,:,:,:), siz(2:end));
    p1     = abs(reshape(input(j,cfg.powindx(:,1),:,:,:), siz(2:end)));
    p2     = abs(reshape(input(j,cfg.powindx(:,2),:,:,:), siz(2:end)));
    
    p      = ipermute(phaseslope(permute(c./sqrt(p1.*p2), pvec), cfg.nbin, cfg.normalize), pvec);
    
    outsum = outsum + p;
    outssq = outssq + p.^2;
  end
  progress('close');
  
elseif length(strfind(cfg.dimord, 'chan'))==2 || length(strfind(cfg.dimord, 'pos'))==2,
  %crossterms are described by chan_chan_therest
  
  siz = size(input);
  if ~hasrpt,
    siz   = [1 siz];
    input = reshape(input, siz);
  end
  
  outsum = zeros(siz(2:end));
  outssq = zeros(siz(2:end));
  pvec   = [3 setdiff(1:numel(siz),3)];
  
  progress('init', cfg.feedback, 'computing metric...');
  for j = 1:siz(1)
    progress(j/siz(1), 'computing metric for replicate %d from %d\n', j, siz(1));
    p1  = zeros([siz(2) 1 siz(4:end)]);
    p2  = zeros([1 siz(3) siz(4:end)]);
    for k = 1:siz(2)
      p1(k,1,:,:,:,:) = input(j,k,k,:,:,:,:);
      p2(1,k,:,:,:,:) = input(j,k,k,:,:,:,:);
    end
    c      = reshape(input(j,:,:,:,:,:,:), siz(2:end));
    p1     = p1(:,ones(1,siz(3)),:,:,:,:);
    p2     = p2(ones(1,siz(2)),:,:,:,:,:);
    p      = ipermute(phaseslope(permute(c./sqrt(p1.*p2),pvec),cfg.nbin, cfg.normalize),pvec);
    outsum = outsum + p;
    outssq = outssq + p.^2;
  end
  progress('close');
  
end

n = siz(1);
c = outsum./n;

if hasrpt,
  if hasjack
    bias = (n-1).^2;
  else
    bias = 1;
  end
  
  v = bias*(outssq - (outsum.^2)./n)./(n - 1);
else
  v = [];
end

%---------------------------------------
function [y] = phaseslope(x, n, norm)

m   = size(x, 1); %total number of frequency bins
y   = zeros(size(x));
x(1:end-1,:,:,:,:) = conj(x(1:end-1,:,:,:,:)).*x(2:end,:,:,:,:);

if strcmp(norm, 'yes')
  coh = zeros(size(x));
  coh(1:end-1,:,:,:,:) = (abs(x(1:end-1,:,:,:,:)) .* abs(x(2:end,:,:,:,:))) + 1;
  %FIXME why the +1? get the coherence
  for k = 1:m
    begindx = max(1,k-n);
    endindx = min(m,k+n);
    y(k,:,:,:,:) = imag(sum(x(begindx:endindx,:,:,:,:)./coh(begindx:endindx,:,:,:,:)));
  end
else
  for k = 1:m
    begindx = max(1,k-n);
    endindx = min(m,k+n);
    y(k,:,:,:,:) = imag(sum(x(begindx:endindx,:,:,:,:)));
  end
end
