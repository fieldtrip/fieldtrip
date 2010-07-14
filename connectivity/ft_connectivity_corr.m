function [c, v, n] = ft_connectivity_corr(input, varargin)

% FIXME write documentation
% takes in a square csd or cov matrix, calculates the partialised
% csd, or partialised cov, if specified, then either normalizes by power (in
% the case of coherence) or returns the partialized csd (if the goal is to
% calculate partial granger, for example.

%FIXME hasrpt should go, because first dim should always be rpt, which could be
%either leave one out or single observations
hasrpt   = keyval('hasrpt',   varargin{:}); if isempty(hasrpt),   hasrpt   = 0;      end
hasjck   = keyval('hasrpt',   varargin{:}); if isempty(hasrpt),   hasrpt   = 0;      end
cmplx    = keyval('complex',  varargin{:}); if isempty(cmplx),    cmplx    = 'abs';  end
feedback = keyval('feedback', varargin{:}); if isempty(feedback), feedback = 'none'; end
dimord   = keyval('dimord',   varargin{:});
powindx  = keyval('powindx',  varargin{:});
pownorm  = keyval('pownorm',  varargin{:}); if isempty(pownorm),  pownorm  = 0;      end
pchanindx   = keyval('pchanindx',   varargin{:});
allchanindx = keyval('allchanindx', varargin{:});

if isempty(dimord)
  error('input parameters should contain a dimord'); 
end

siz = size(input);
if ~hasrpt,
  siz   = [1 siz];
  input = reshape(input, siz);
end

% do partialisation if necessary
if ~isempty(pchanindx),
  % partial spectra are computed as in Rosenberg JR et al (1998) J.
  % Neuroscience Methods, equation 38
  
  chan   = allchanindx;
  nchan  = numel(chan);
  pchan  = pchanindx;
  npchan = numel(pchan);
  newsiz = siz;
  newsiz(2:3) = numel(chan); % size of partialised csd
  
  A  = zeros(newsiz);
  
  % FIXME this only works for data without time dimension
  if numel(siz)>4, error('this only works for data without time'); end
  for j = 1:siz(1) %rpt loop
    AA = reshape(input(j, chan,  chan, : ), [nchan  nchan  siz(4:end)]);
    AB = reshape(input(j, chan,  pchan,: ), [nchan  npchan siz(4:end)]);
    BA = reshape(input(j, pchan, chan, : ), [npchan nchan  siz(4:end)]);
    BB = reshape(input(j, pchan, pchan, :), [npchan npchan siz(4:end)]);
    for k = 1:siz(4) %freq loop
      A(j,:,:,k) = AA(:,:,k) - AB(:,:,k)*pinv(BB(:,:,k))*BA(:,:,k);
    end
  end
  input = A;
  siz = size(input);
else
  % do nothing
end

if (length(strfind(dimord, 'chan'))~=2 || length(strfind(dimord, 'pos'))>0) && ~isempty(powindx),
  %crossterms are not described with chan_chan_therest, but are linearly indexed
  outsum = zeros(siz(2:end));
  outssq = zeros(siz(2:end));
  
  progress('init', feedback, 'computing metric...');
  for j = 1:siz(1)
    progress(j/siz(1), 'computing metric for replicate %d from %d\n', j, siz(1));
    if pownorm
      p1    = reshape(input(j,powindx(:,1),:,:,:), siz(2:end));
      p2    = reshape(input(j,powindx(:,2),:,:,:), siz(2:end));
      denom = sqrt(p1.*p2); clear p1 p2
    else
      denom = 1;
    end
    outsum = outsum + complexeval(reshape(input(j,:,:,:,:), siz(2:end))./denom, cmplx);
    outssq = outssq + complexeval(reshape(input(j,:,:,:,:), siz(2:end))./denom, cmplx).^2;
  end
  progress('close');
  
elseif length(strfind(dimord, 'chan'))==2 || length(strfind(dimord, 'pos'))==2,
  %crossterms are described by chan_chan_therest
  
  outsum = zeros(siz(2:end));
  outssq = zeros(siz(2:end));
  progress('init', feedback, 'computing metric...');
  for j = 1:siz(1)
    progress(j/siz(1), 'computing metric for replicate %d from %d\n', j, siz(1));
    if pownorm
      p1  = zeros([siz(2) 1 siz(4:end)]);
      p2  = zeros([1 siz(3) siz(4:end)]);
      for k = 1:siz(2)
        p1(k,1,:,:,:,:) = input(j,k,k,:,:,:,:);
        p2(1,k,:,:,:,:) = input(j,k,k,:,:,:,:);
      end
      p1    = p1(:,ones(1,siz(3)),:,:,:,:);
      p2    = p2(ones(1,siz(2)),:,:,:,:,:);
      denom = sqrt(p1.*p2); clear p1 p2;
    else
      denom = 1;
    end
    outsum = outsum + complexeval(reshape(input(j,:,:,:,:,:,:), siz(2:end))./denom, cmplx);
    outssq = outssq + complexeval(reshape(input(j,:,:,:,:,:,:), siz(2:end))./denom, cmplx).^2;
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

function [c] = complexeval(c, str)

switch str
  case 'complex'
    %do nothing
  case 'abs'
    c = abs(c);
  case 'angle'
    c = angle(c);
  case 'imag'
    c = imag(c);
  case 'real'
    c = real(c);
  otherwise
    error('complex = ''%s'' not supported', cmplx);
end
