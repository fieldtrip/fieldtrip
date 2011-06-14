function [granger, v, n] = ft_connectivity_granger(H, Z, S, fs, hasjack, powindx)

%Usage: causality = hz2causality(H,S,Z,fs);
%Inputs: transfer  = transfer function,
%        crsspctrm = 3-D spectral matrix;
%        noisecov  = noise covariance,
%        fs        = sampling rate
%Outputs: granger (Granger causality between all channels)
%               : auto-causality spectra are set to zero
% Reference: Brovelli, et. al., PNAS 101, 9849-9854 (2004).
%M. Dhamala, UF, August 2006.

%FIXME speed up code and check
siz = size(H);
if numel(siz)==4,
  siz(5) = 1;
end
n   = siz(1);
Nc  = siz(2);

outsum = zeros(siz(2:end));
outssq = zeros(siz(2:end));

if isempty(powindx),
  % data are chan_chan_therest
  for kk = 1:n
    for ii = 1:Nc
      for jj = 1:Nc
        if ii ~=jj,
          zc     = reshape(Z(kk,jj,jj,:) - Z(kk,ii,jj,:).^2./Z(kk,ii,ii,:),[1 1 1 1 siz(5)]);
          zc     = repmat(zc,[1 1 1 siz(4) 1]);
          numer  = reshape(abs(S(kk,ii,ii,:,:)),[1 1 siz(4:end)]);
          denom  = reshape(abs(S(kk,ii,ii,:,:)-zc.*abs(H(kk,ii,jj,:,:)).^2./fs),[1 1 siz(4:end)]);
          outsum(jj,ii,:,:) = outsum(jj,ii,:,:) + log(numer./denom);
          outssq(jj,ii,:,:) = outssq(jj,ii,:,:) + (log(numer./denom)).^2;
        end
      end
      outsum(ii,ii,:,:) = 0;%self-granger set to zero
    end
  end
elseif ~iscell(powindx) && ~isstruct(powindx)
  % data are linearly indexed
  for j = 1:n
    for k = 1:Nc
      %iauto1  = sum(powindx==powindx(k,1),2)==2;
      %iauto2  = sum(powindx==powindx(k,2),2)==2;
      %icross1 = k;
      %icross2 = sum(powindx==powindx(ones(Nc,1)*k,[2 1]),2)==2;
      if mod(k-1, 4)==0
        iauto1=k;iauto2=k;icross1=k;icross2=k;
      elseif mod(k-1, 4)==1
        iauto1=k+2;iauto2=k-1;icross1=k;icross2=k+1;
      elseif mod(k-1, 4)==2
        iauto1=k-2;iauto2=k+1;icross1=k;icross2=k-1;
      elseif mod(k-1, 4)==3
        iauto1=k;iauto2=k;icross1=k;icross2=k;
      end
      
      zc      = Z(j,iauto2,:,:) - Z(j,icross1,:,:).^2./Z(j,iauto1,:,:);
      numer   = abs(S(j,iauto1,:,:));
      denom   = abs(S(j,iauto1,:,:)-zc(:,:,ones(1,size(H,3)),:).*abs(H(j,icross1,:,:)).^2./fs);
      outsum(icross2,:,:) = outsum(icross2,:,:) + reshape(log(numer./denom), [1 siz(3:end)]);
      outssq(icross2,:,:) = outssq(icross2,:,:) + reshape((log(numer./denom)).^2, [1 siz(3:end)]);
    end
  end
elseif iscell(powindx)
  % blockwise granger
  % H = transfer function nchan x nchan x nfreq
  % Z = noise covariance  nchan x nchan
  % S = crosspectrum      nchan x nchan x nfreq
  % powindx{1} is a list of indices for block1
  % powindx{2} is a list of indices for block2
  
  %FIXME rewrite to allow for multiple blocks
  %FIXME change cfg.block functionality in this case
  %cfg.blockindx = {{list of channel names} [list of block indices]}
  block1 = powindx{1}(:);
  block2 = powindx{2}(:);
  
  n     = size(H,1);
  nchan = size(H,2);
  nfreq = size(H,4);
  
  n1 = numel(block1);
  n2 = numel(block2);
  
  % reorder
  S = S(:,[block1;block2],[block1;block2],:);
  H = H(:,[block1;block2],[block1;block2],:);
  Z = Z(:,[block1;block2],[block1;block2]);
  
  indx1 = 1:n1;
  indx2 = (n1+1):(n1+n2);
  
  outsum = zeros(2,2,nfreq);
  outssq = zeros(2,2,nfreq);
  for k = 1:n
    tmpZ = reshape(Z(k,:,:), [nchan nchan]);
    
    % projection matrix for block2 -> block1
    P1 = [eye(n1)                                zeros(n1,n2);
      -tmpZ(indx2,indx1)/tmpZ(indx1,indx1)     eye(n2)];
    
    % projection matrix for block1 -> block2
    P2 = [  eye(n1)    -tmpZ(indx1,indx2)/tmpZ(indx2,indx2);
      zeros(n2,n1) eye(n2)];
    
    % invert only once
    invP1 = inv(P1);
    invP2 = inv(P2);
    for jj = 1:nfreq
      % post multiply transfer matrix with the inverse of the projection matrix
      % this is equivalent to time domain pre multiplication with P
      Sj = reshape(S(k,:,:,jj), [nchan nchan]);
      Zj = tmpZ(:,:);
      H1 = reshape(H(k,:,:,jj), [nchan nchan])*invP1;
      H2 = reshape(H(k,:,:,jj), [nchan nchan])*invP2;
      num1 = abs(det(Sj(indx1,indx1))); % numerical round off leads to tiny imaginary components
      num2 = abs(det(Sj(indx2,indx2))); % numerical round off leads to tiny imaginary components
      denom1 = abs(det(H1(indx1,indx1)*Zj(indx1,indx1)*H1(indx1,indx1)'));
      denom2 = abs(det(H2(indx2,indx2)*Zj(indx2,indx2)*H2(indx2,indx2)'));
      %rH1 = real(H1(indx1,indx1));
      %rH2 = real(H2(indx2,indx2));
      %iH1 = imag(H1(indx1,indx1));
      %iH2 = imag(H2(indx2,indx2));
      %h1 = rH1*Zj(indx1,indx1)*rH1' + iH1*Zj(indx1,indx1)*iH1';
      %h2 = rH2*Zj(indx2,indx2)*rH2' + iH2*Zj(indx2,indx2)*iH2';
      %denom1 = abs(det(h1));
      %denom2 = abs(det(h2));
      
      outsum(2,1,jj) = log( num1./denom1 )    + outsum(2,1,jj);
      outsum(1,2,jj) = log( num2./denom2 )    + outsum(1,2,jj);
      outssq(2,1,jj) = log( num1./denom1 ).^2 + outssq(2,1,jj);
      outssq(1,2,jj) = log( num2./denom2 ).^2 + outssq(1,2,jj);
    end
  end
elseif isstruct(powindx)
  %blockwise conditional
  
  n     = size(H,1);
  ncmb  = size(H,2);
  nfreq = size(H,3);
  ncnd  = size(powindx.cmbindx,1);
  
  outsum = zeros(ncnd, nfreq);
  outssq = zeros(ncnd, nfreq);
  for k = 1:n
    tmpS = reshape(S, [ncmb nfreq]);
    tmpH = reshape(H, [ncmb nfreq]);
    tmpZ = reshape(Z, [ncmb 1]);
    tmp  = blockwise_conditionalgranger(tmpS,tmpH,tmpZ,powindx.cmbindx,powindx.n);
    
    outsum = outsum + tmp;
    outssq = outssq + tmp.^2;
  end
end

granger = outsum./n;
if n>1,
  if hasjack
    bias = (n-1).^2;
  else
    bias = 1;
  end
  v = bias*(outssq - (outsum.^2)./n)./(n - 1);
else
  v = [];
end
