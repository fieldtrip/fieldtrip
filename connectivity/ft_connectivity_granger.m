function [granger, v, n] = ft_connectivity_granger(H, Z, S, varargin)

% FT_CONNECTIVITY_GRANGER computes spectrally resolved granger causality. This
% implementation is loosely based on the code used in Brovelli, et. al., PNAS 101,
% 9849-9854 (2004).
%
% Use as
%   [GRANGER, V, N] = FT_CONNECTIVITY_GRANGER(H, Z, S, ...)
%
% The input data should be
%   H = spectral transfer matrix, Nrpt x Nchan x Nchan x Nfreq (x Ntime),
%       or Nrpt x Nchancmb x Nfreq (x Ntime). Nrpt can be 1.
%   Z = the covariance matrix of the noise, Nrpt x Nchan x Nchan (x Ntime),
%       or Nrpt x Nchancmb (x Ntime).
%   S = the cross-spectral density matrix with the same dimensionality as H.
%
% Additional optional input arguments come as key-value pairs:
%   'dimord'  = required string specifying how to interpret the input data
%               supported values are 'rpt_chan_chan_freq(_time) and
%               'rpt_chan_freq(_time), 'rpt_pos_pos_XXX' and 'rpt_pos_XXX'
%   'method'  = 'granger' (default), or 'instantaneous', or 'total'.
%   'hasjack' = 0 (default) is a boolean specifying whether the input
%               contains leave-one-outs, required for correct variance
%               estimate
%   'powindx' = is a variable determining the exact computation, see below
%
% If the inputdata is such that the channel-pairs are linearly indexed, granger
% causality is computed per quadruplet of consecutive entries, where the convention
% is as follows:
%
%  H(:, (k-1)*4 + 1, :, :, :) -> 'chan1-chan1'
%  H(:, (k-1)*4 + 2, :, :, :) -> 'chan1->chan2'
%  H(:, (k-1)*4 + 3, :, :, :) -> 'chan2->chan1'
%  H(:, (k-1)*4 + 4, :, :, :) -> 'chan2->chan2'
%
% The same holds for the Z and S matrices.
%
% Pairwise block-granger causality can be computed when the inputdata has
% dimensionality Nchan x Nchan. In that case powindx should be specified, as a 1x2
% cell-array indexing the individual channels that go into each 'block'.
%
% See also FT_CONNECTIVITYANALYSIS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Undocumented option: powindx can be a struct. In that case, blockwise
% conditional granger can be computed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2009-2017, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

method  = ft_getopt(varargin, 'method',  'granger');
hasjack = ft_getopt(varargin, 'hasjack', 0);
powindx = ft_getopt(varargin, 'powindx', []);
dimord  = ft_getopt(varargin, 'dimord',  []);

%FIXME speed up code and check
siz = size(H);
if numel(siz)==4
  siz(5) = 1;
end
n   = siz(1);
Nc  = siz(2);

outsum = zeros(siz(2:end));
outssq = zeros(siz(2:end));

% crossterms are described by chan_chan_therest
issquare = length(strfind(dimord, 'chan'))==2 || length(strfind(dimord, 'pos'))==2;

switch method
  case 'granger'
    
    if issquare && isempty(powindx)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % data are chan_chan_therest
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      for kk = 1:n
        for ii = 1:Nc
          for jj = 1:Nc
            if ii ~=jj
              zc     = reshape(Z(kk,jj,jj,:) - Z(kk,ii,jj,:).^2./Z(kk,ii,ii,:),[1 1 1 1 siz(5)]);
              zc     = repmat(zc,[1 1 1 siz(4) 1]);
              numer  = reshape(abs(S(kk,ii,ii,:,:)),[1 1 siz(4:end)]);
              denom  = reshape(abs(S(kk,ii,ii,:,:)-zc.*abs(H(kk,ii,jj,:,:)).^2),[1 1 siz(4:end)]);
              outsum(jj,ii,:,:) = outsum(jj,ii,:,:) + log(numer./denom);
              outssq(jj,ii,:,:) = outssq(jj,ii,:,:) + (log(numer./denom)).^2;
            end
          end
          outsum(ii,ii,:,:) = 0;%self-granger set to zero
        end
      end
      
    elseif ~issquare && isempty(powindx)
      %%%%%%%%%%%%%%%%%%%%%%%%%%
      %data are linearly indexed
      %%%%%%%%%%%%%%%%%%%%%%%%%%
      
      for j = 1:n
        for k = 1:Nc
          %FIXME powindx is not used here anymore
          %iauto1  = sum(powindx==powindx(k,1),2)==2;
          %iauto2  = sum(powindx==powindx(k,2),2)==2;
          %icross1 = k;
          %icross2 = sum(powindx==powindx(ones(Nc,1)*k,[2 1]),2)==2;
          
          % The following is based on hard-coded assumptions: which is fair
          % to do if the order of the labelcmb is according to the output of
          % ft_connectivity_csd2transfer
          if mod(k-1, 4)==0
            continue; % auto granger set to 0
            %iauto1=k;iauto2=k;icross1=k;icross2=k;
          elseif mod(k-1, 4)==1
            iauto1=k+2;iauto2=k-1;icross1=k;icross2=k+1;
          elseif mod(k-1, 4)==2
            iauto1=k-2;iauto2=k+1;icross1=k;icross2=k-1;
          elseif mod(k-1, 4)==3
            continue; % auto granger set to 0
            %iauto1=k;iauto2=k;icross1=k;icross2=k;
          end
          
          zc      = Z(j,iauto2,:,:) - Z(j,icross1,:,:).^2./Z(j,iauto1,:,:);
          numer   = abs(S(j,iauto1,:,:));
          denom   = abs(S(j,iauto1,:,:)-zc(:,:,ones(1,size(H,3)),:).*abs(H(j,icross1,:,:)).^2);
          outsum(icross2,:,:) = outsum(icross2,:,:) + reshape(log(numer./denom), [1 siz(3:end)]);
          outssq(icross2,:,:) = outssq(icross2,:,:) + reshape((log(numer./denom)).^2, [1 siz(3:end)]);
        end
      end
      
    elseif issquare && iscell(powindx)
      %%%%%%%%%%%%%%%%%%%
      % blockwise granger
      %%%%%%%%%%%%%%%%%%%
      
      % H = transfer function nchan x nchan x nfreq
      % Z = noise covariance  nchan x nchan
      % S = crosspectrum      nchan x nchan x nfreq
      % powindx{k} is a list of indices for block k
      
      nblock = numel(powindx);
      n      = size(H,1);
      nfreq  = size(H,4);
      
      outsum = zeros(nblock,nblock,nfreq);
      outssq = zeros(nblock,nblock,nfreq);
      
      for k = 1:nblock
        for m = (k+1):nblock
          indx  = [powindx{k}(:);powindx{m}(:)];
          indx  = cat(1,indx,setdiff((1:size(Z,2))',indx));
          n1    = numel(powindx{k});
          n2    = numel(powindx{m});
          ntot  = numel(indx);
          indx1 = 1:n1;
          indx1_  = (1:ntot)'; indx1_(indx1) = [];
          indx2   = n1+(1:n2);
          indx12_ = (1:ntot)'; indx12_([indx1(:);indx2(:)]) = []; 
          
          for kk = 1:n
            tmpZ = reshape(Z(kk,indx,indx), [ntot ntot]);
            
            % projection matrix for therest+block2 -> block1
            P1 = [eye(n1)                                zeros(n1,ntot-n1);
              -tmpZ(indx1_,indx1)/tmpZ(indx1,indx1)      eye(ntot-n1)];
            
            % projection matrix for therest+block1 -> block2
            P2 = [  eye(n1)    -tmpZ(indx1,indx2)/tmpZ(indx2,indx2) zeros(n1,ntot-n1-n2);
                  zeros(n2,n1)   eye(n2) zeros(n2, ntot-n1-n2);
                  zeros(ntot-n1-n2,n1) -tmpZ(indx12_,indx2)/tmpZ(indx2,indx2) eye(ntot-n1-n2)];
            
            % invert only once
            for jj = 1:nfreq
              % post multiply transfer matrix with the inverse of the projection matrix
              % this is equivalent to time domain pre multiplication with P
              Sj = reshape(S(kk,indx,indx,jj), [ntot ntot]);
              Zj = tmpZ;%(:,:);
              H1 = reshape(H(kk,indx,indx,jj), [ntot ntot])/P1;
              H2 = reshape(H(kk,indx,indx,jj), [ntot ntot])/P2;
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
              
              outsum(m,k,jj) = log( num1./denom1 )    + outsum(m,k,jj);
              outsum(k,m,jj) = log( num2./denom2 )    + outsum(k,m,jj);
              outssq(m,k,jj) = log( num1./denom1 ).^2 + outssq(m,k,jj);
              outssq(k,m,jj) = log( num2./denom2 ).^2 + outssq(k,m,jj);
            end
          end
          
        end
      end
      
    elseif ~issquare && isstruct(powindx) && isfield(powindx, 'n')
      %%%%%%%%%%%%%%%%%%%%%%
      %blockwise conditional
      %%%%%%%%%%%%%%%%%%%%%%
      
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
      
    elseif ~issquare && isstruct(powindx)
      %%%%%%%%%%%%%%%%%%%%%%
      %triplet conditional
      %%%%%%%%%%%%%%%%%%%%%%
      
      % decode from the powindx struct which rows in the data correspond with
      % the triplets, and which correspond with the duplets
      ublockindx = unique(powindx.blockindx);
      nperblock  = zeros(size(ublockindx));
      for k = 1:numel(ublockindx)
        nperblock(k,1) = sum(powindx.blockindx==ublockindx(k));
      end
      if ~all(ismember(nperblock,[4 9]))
        error('the data should be a mixture of trivariate and bivariate decompositions');
      end
      indx_triplets = ismember(powindx.blockindx, ublockindx(nperblock==9)); ntriplets = sum(indx_triplets)./9;
      indx_duplets  = ismember(powindx.blockindx, ublockindx(nperblock==4)); nduplets  = sum(indx_duplets)./4;
      
      % this assumes well-behaved powindx.cmbindx
      cmbindx2 = reshape(powindx.cmbindx(indx_duplets, 1), 2, [])';
      cmbindx3 = reshape(powindx.cmbindx(indx_triplets,1), 3, [])';
      
      cmbindx2 = cmbindx2(1:2:end,:);
      cmbindx3 = cmbindx3(1:3:end,:);
      
      cmbindx  = powindx.outindx;
      
      n   = size(H,1);
      siz = size(H);
      
      outsum = zeros(size(cmbindx,1), size(H,3));
      outssq = outsum;
      
      % call the low-level function
      for k = 1:n
        H3 = reshape(H(k,indx_triplets,:,:), [3 3 ntriplets siz(3:end)]);
        Z3 = reshape(Z(k,indx_triplets),     [3 3 ntriplets]);
        
        H2 = reshape(H(k,indx_duplets,:,:), [2 2 nduplets siz(3:end)]);
        Z2 = reshape(Z(k,indx_duplets),     [2 2 nduplets]);
        
        tmp  = triplet_conditionalgranger(H3,Z3,cmbindx3,H2,Z2,cmbindx2,cmbindx);
        
        outsum = outsum + tmp;
        outssq = outssq + tmp.^2;
      end
      
      
      
    end
    
  case 'instantaneous'
    
    if issquare && isempty(powindx)
      % data are chan_chan_therest
      for kk = 1:n
        for ii = 1:Nc
          for jj = 1:Nc
            if ii ~=jj
              zc1    = reshape(Z(kk,jj,jj,:) - Z(kk,ii,jj,:).^2./Z(kk,ii,ii,:),[1 1 1 1 siz(5)]);
              zc2    = reshape(Z(kk,ii,ii,:) - Z(kk,jj,ii,:).^2./Z(kk,jj,jj,:),[1 1 1 1 siz(5)]);
              zc1    = repmat(zc1,[1 1 1 siz(4) 1]);
              zc2    = repmat(zc2,[1 1 1 siz(4) 1]);
              term1  = abs(S(kk,ii,ii,:,:)) - zc1.*abs(H(kk,ii,jj,:,:)).^2;
              term2  = abs(S(kk,jj,jj,:,:)) - zc2.*abs(H(kk,jj,ii,:,:)).^2;
              numer  = term1.*term2;
              denom  = abs(S(kk,ii,ii,:,:).*S(kk,jj,jj,:,:) - S(kk,ii,jj,:,:).*S(kk,jj,ii,:,:));
              outsum(jj,ii,:,:) = outsum(jj,ii,:,:) + reshape(log(numer./denom),  [1 1 siz(4:end)]);
              outssq(jj,ii,:,:) = outssq(jj,ii,:,:) + reshape((log(numer./denom)).^2, [1 1 siz(4:end)]);
            end
          end
          outsum(ii,ii,:,:) = 0;%self-granger set to zero
        end
      end
    elseif ~issquare && isempty(powindx)
      % data are linearly indexed
      for j = 1:n
        for k = 1:Nc
          %iauto1  = sum(powindx==powindx(k,1),2)==2;
          %iauto2  = sum(powindx==powindx(k,2),2)==2;
          %icross1 = k;
          %icross2 = sum(powindx==powindx(ones(Nc,1)*k,[2 1]),2)==2;
          if mod(k-1, 4)==0
            continue; % auto granger set to 0
            %iauto1=k;iauto2=k;icross1=k;icross2=k;
          elseif mod(k-1, 4)==1
            iauto1=k+2;iauto2=k-1;icross1=k;icross2=k+1;
          elseif mod(k-1, 4)==2
            iauto1=k-2;iauto2=k+1;icross1=k;icross2=k-1;
          elseif mod(k-1, 4)==3
            continue; % auto granger set to 0
            %iauto1=k;iauto2=k;icross1=k;icross2=k;
          end
          
          zc1     = Z(j,iauto1,:, :) - Z(j,icross2,:, :).^2./Z(j,iauto2,:, :);
          zc2     = Z(j,iauto2,:, :) - Z(j,icross1,:, :).^2./Z(j,iauto1,:, :);
          term1   = abs(S(j,iauto2,:,:)) - zc1(:,:,ones(1,size(H,3)),:).*abs(H(j,icross2,:,:)).^2;
          term2   = abs(S(j,iauto1,:,:)) - zc2(:,:,ones(1,size(H,3)),:).*abs(H(j,icross1,:,:)).^2;
          numer   = term1.*term2;
          denom   = abs(S(j,iauto1,:,:).*S(j,iauto2,:,:) - S(j,icross1,:,:).*S(j,icross2,:,:));
          
          outsum(icross2,:,:) = outsum(icross2,:,:) + reshape(log(numer./denom), [1 siz(3:end)]);
          outssq(icross2,:,:) = outssq(icross2,:,:) + reshape((log(numer./denom)).^2, [1 siz(3:end)]);
        end
      end
    elseif issquare && iscell(powindx)
      % blockwise granger
      % H = transfer function nchan x nchan x nfreq
      % Z = noise covariance  nchan x nchan
      % S = crosspectrum      nchan x nchan x nfreq
      % powindx{1} is a list of indices for block1
      % powindx{2} is a list of indices for block2
      ft_error('instantaneous causality is not implemented for blockwise factorizations');
    elseif isstruct(powindx)
      %blockwise conditional
      ft_error('blockwise conditional instantaneous causality is not implemented');
    else
      ft_error('not implemented');
    end
    
  case 'total'
    
    if issquare && isempty(powindx)
      % data are chan_chan_therest
      for kk = 1:n
        for ii = 1:Nc
          for jj = 1:Nc
            if ii ~=jj
              numer  = abs(S(kk,ii,ii,:,:).*S(kk,jj,jj,:,:));
              denom  = abs(S(kk,ii,ii,:,:).*S(kk,jj,jj,:,:) - S(kk,ii,jj,:,:).*S(kk,jj,ii,:,:));
              outsum(jj,ii,:,:) = outsum(jj,ii,:,:) + reshape(log(numer./denom), [1 1 siz(4:end)]);
              outssq(jj,ii,:,:) = outssq(jj,ii,:,:) + reshape((log(numer./denom)).^2, [1 1 siz(4:end)]);
            end
          end
          outsum(ii,ii,:,:) = 0;%self-granger set to zero
        end
      end
    elseif ~issquare && isempty(powindx)
      % data are linearly indexed
      for j = 1:n
        for k = 1:Nc
          %iauto1  = sum(powindx==powindx(k,1),2)==2;
          %iauto2  = sum(powindx==powindx(k,2),2)==2;
          %icross1 = k;
          %icross2 = sum(powindx==powindx(ones(Nc,1)*k,[2 1]),2)==2;
          if mod(k-1, 4)==0
            continue; % auto granger set to 0
            %iauto1=k;iauto2=k;icross1=k;icross2=k;
          elseif mod(k-1, 4)==1
            iauto1=k+2;iauto2=k-1;icross1=k;icross2=k+1;
          elseif mod(k-1, 4)==2
            iauto1=k-2;iauto2=k+1;icross1=k;icross2=k-1;
          elseif mod(k-1, 4)==3
            continue; % auto granger set to 0
            %iauto1=k;iauto2=k;icross1=k;icross2=k;
          end
          
          numer   = abs(S(j,iauto1,:,:).*S(j,iauto2,:,:));
          denom   = abs(S(j,iauto1,:,:).*S(j,iauto2,:,:) - S(j,icross1,:,:).*S(j,icross2,:,:));
          outsum(icross2,:,:) = outsum(icross2,:,:) + reshape(log(numer./denom), [1 siz(3:end)]);
          outssq(icross2,:,:) = outssq(icross2,:,:) + reshape((log(numer./denom)).^2, [1 siz(3:end)]);
        end
      end
    elseif issquare && iscell(powindx)
      % blockwise total interdependence
      
      % S = crosspectrum      nchan x nchan x nfreq
      % powindx{k} is a list of indices for block k
      
      nblock = numel(powindx);
      n      = size(S,1);
      nfreq  = size(S,4);
      
      outsum = zeros(nblock,nblock,nfreq);
      outssq = zeros(nblock,nblock,nfreq);
      
      for k = 1:nblock
        for m = (k+1):nblock
          indx  = [powindx{k}(:);powindx{m}(:)];
          n1    = numel(powindx{k});
          n2    = numel(powindx{m});
          ntot  = n1+n2;
          indx1 = 1:n1;
          indx2 = (n1+1):ntot;
          
          for kk = 1:n
            for jj = 1:nfreq
              Sj = reshape(S(kk,indx,indx,jj), [ntot ntot]);
              num1 = abs(det(Sj(indx1,indx1))); % numerical round off leads to tiny imaginary components
              num2 = abs(det(Sj(indx2,indx2))); % numerical round off leads to tiny imaginary components
              num  = num1.*num2;
              denom = abs(det(Sj));
              
              outsum(m,k,jj) = log( num./denom )    + outsum(m,k,jj);
              outsum(k,m,jj) = log( num./denom )    + outsum(k,m,jj);
              outssq(m,k,jj) = log( num./denom ).^2 + outssq(m,k,jj);
              outssq(k,m,jj) = log( num./denom ).^2 + outssq(k,m,jj);
            end
          end
          
        end
      end
      
    elseif issquare && isstruct(powindx)
      %blockwise conditional
      ft_error('blockwise conditional total interdependence is not implemented');
    end
    
  case 'iis'
    ft_warning('THIS IS EXPERIMENTAL CODE, USE AT YOUR OWN RISK!');
    % this is experimental
    if ~issquare && isempty(powindx)
      A = transfer2coeffs(shiftdim(H),(0:size(H,3)-1));
      ncmb = size(A,1)./4;
      iis = coeffs2iis(reshape(A,[2 2 ncmb size(A,2)]),reshape(Z,[2 2 ncmb]));
      iis = repmat(iis(:),[1 4])';
      outsum = iis(:);
      outssq = nan(size(outsum));
    else
      ft_error('iis can only be computed when the input contains sets of bivariate factorizations');
    end
    
  otherwise
    ft_error('unsupported output requested');
end

granger = outsum./n;
if n>1
  if hasjack
    bias = (n-1).^2;
  else
    bias = 1;
  end
  v = bias*(outssq - (outsum.^2)./n)./(n - 1);
else
  v = [];
end
