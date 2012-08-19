function [output] = ft_connectivity_csd2transfer(freq, varargin)

% CSD2TRANSFER computes the transfer-function from frequency domain data
% using the Wilson-Burg algorithm. The transfer function can be used for
% the computation of directional measures of connectivity, such as granger
% causality, partial directed coherence, or directed transfer functions
%
% Use as
%   [output] = ft_connectivity_csd2transfer(freq, varargin)
%
% Where freq is a data structure containing frequency domain data containing
% the cross-spectral density computed between all pairs of channels, thus
% containing a 'dimord' of 'chan_chan_freq(_time)'.
%
% Configurable options come in key-value pairs:
%
%   numiteration = scalar value (default: 100) the number of iterations
%   channelcmb   = Nx2 cell-array listing the channel pairs for the spectral
%                    factorization. If not defined or empty (default), a
%                    full multivariate factorization is performed, otherwise
%                    a multiple pairwise factorization is done.
%   tol          = scalar value (default: 1e-18) tolerance limit truncating
%                    the iterations
%
% The code for the Wilson-Burg algorithm has been very generously provided by 
% Dr. Mukesh Dhamala, and Prof. Mingzhou Ding and his group. 
%
% If you use this code for studying directed interactions, please cite from
% the following references:
% -M.Dhamala, G.Rangarajan, M.Ding, Physical Review Letters 100, 018701 (2008)
% -M.Dhamala, G.rangarajan, M.Ding, Neuroimage 41, 354 (2008)

% Undocumented options:
%
%   block
%   blockindx
%   svd
%
% Copyright (C) 2009-2011, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

numiteration = ft_getopt(varargin, 'numiteration', 100);
channelcmb   = ft_getopt(varargin, 'channelcmb',   {});
block        = ft_getopt(varargin, 'block',        {});
blockindx    = ft_getopt(varargin, 'blockindx',    cell(0,2));
tol          = ft_getopt(varargin, 'tol',          1e-18);
fb           = ft_getopt(varargin, 'feedback',     'textbar');
sfmethod     = ft_getopt(varargin, 'sfmethod',     'multivariate');
dosvd        = ft_getopt(varargin, 'svd',          'no');

dosvd = istrue(dosvd);

% check whether input data is valid
freq = ft_checkdata(freq, 'datatype', 'freq');
if ~isfield(freq, 'crsspctrm') || ~isfield(freq, 'label')
  error('the input data does not contain cross-spectral density data in the supported format');
end

hasrpt = ~isempty(strfind(freq.dimord, 'rpt'));
if hasrpt,
  nrpt = numel(freq.cumtapcnt);
else
  nrpt = 1;
end

if ~isempty(channelcmb) && ~(strcmp(channelcmb{1}, 'all') && strcmp(channelcmb{2}, 'all')) && ~strcmp(sfmethod, 'bivariate')
  % enforce sfmethod to be bivariate
  sfmethod = 'bivariate';
end

% if bivariate is requested without channelcmb, do all versus all pairwise
if isempty(channelcmb) && strcmp(sfmethod, 'bivariate')
  channelcmb = {'all' 'all'};
end


if ~isempty(channelcmb)
  if numel(channelcmb)==2 && strcmp(channelcmb{1},'all') && strcmp(channelcmb{2}, 'all')
    [cmbindx(:,1), cmbindx(:,2)] = find(tril(ones(numel(freq.label)),-1));
    ok = true(size(cmbindx,1),1);
  end
  channelcmb = ft_channelcombination(channelcmb, freq.label);
end

if ~isempty(block)
  % sanity check 1
  if ~iscell(block) 
    error('block should be a cell-array containing 3 cells');
  end
  % sanity check 2
  if strcmp(sfmethod, 'bivariate')
    error('when block is specified, it is not OK to do bivariate decomposition');
  end
end

if isfield(freq, 'time'),
  ntim = numel(freq.time);
else
  ntim = 1;
end

siz = size(freq.crsspctrm);
if ntim==1,
  siz = [siz 1]; %add dummy dimensionality for time axis
end

if strcmp(sfmethod, 'bivariate')
  list1 = char(channelcmb(:,1)); siz1 = size(list1,2);
  list2 = char(channelcmb(:,2)); siz2 = size(list2,2);
  ix    = 1:min(siz1,siz2);
  npair = sum(sum(list1(:,ix)~=list2(:,ix),2)~=0);
  fprintf('computing pairwise non-parametric spectral factorization on %d channel pairs\n', npair);%size(channelcmb,1) - numel(unique(channelcmb(:))));
elseif strcmp(sfmethod, 'multivariate')
  fprintf('computing multivariate non-parametric spectral factorization on %d channels\n', numel(freq.label));   
else
  error('unknown sfmethod %s', sfmethod);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computations start here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(sfmethod, 'multivariate') && isempty(block) && nrpt>1,
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % multiple repetitions, loop over repetitions
  % multivariate decomposition
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  H = zeros(siz) + 1i.*zeros(siz);
  S = zeros(siz) + 1i.*zeros(siz);
  Z = zeros([siz(1:3) siz(end)]);
  for k = 1:nrpt
    for m = 1:ntim
      tmp = reshape(freq.crsspctrm(k,:,:,:,m), siz(2:end-1));
      [Htmp, Ztmp, Stmp] = sfactorization_wilson(tmp, freq.freq, ...
                                                   numiteration, tol, fb);
      H(k,:,:,:,m) = Htmp;
      Z(k,:,:,m)   = Ztmp;
      S(k,:,:,:,m) = Stmp;
    end 
  end
  
elseif strcmp(sfmethod, 'multivariate') && isempty(block),
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % standard code
  % multivariate decomposition
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  H   = zeros(siz) + 1i.*zeros(siz);
  S   = zeros(siz) + 1i.*zeros(siz);
  Z   = zeros([siz(1:2) siz(end)]);
  
  % only do decomposition once
  for m = 1:ntim
    tmp = freq.crsspctrm(:,:,:,m);

    % do SVD to avoid zigzags due to numerical issues
    if dosvd
      dat     = sum(tmp,3);
      [u,s,v] = svd(real(dat));
      for k = 1:size(tmp,3)
        tmp(:,:,k) = u'*tmp(:,:,k)*u;
      end
    end
    
    if any(isnan(tmp(:))),
      Htmp = nan;
      Ztmp = nan;
      Stmp = nan;
    else
      [Htmp, Ztmp, Stmp] = sfactorization_wilson(tmp, freq.freq, ...
                                                   numiteration, tol, fb);
    end
    
    % undo SVD
    if dosvd
      for k = 1:size(tmp,3)
        Htmp(:,:,k) = u*Htmp(:,:,k)*u';
        Stmp(:,:,k) = u*Stmp(:,:,k)*u';
      end
      Ztmp = u*Ztmp*u';
    end
    
    H(:,:,:,m) = Htmp;
    Z(:,:,m)   = Ztmp;
    S(:,:,:,m) = Stmp;
  end
elseif strcmp(sfmethod, 'bivariate') && nrpt>1,
  % error 
  error('single trial estimates and linear combination indexing is not implemented');
  
elseif nrpt>1 && ~isempty(block),
  % error 
  error('single trial estimates and blockwise factorisation is not yet implemented');
  
elseif strcmp(sfmethod, 'bivariate')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % pairwise factorization resulting in linearly indexed transfer functions
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %convert list of channel labels into indices
  if ~exist('cmbindx', 'var')
    cmbindx     = zeros(size(channelcmb));
    ok          = true(size(cmbindx,1), 1);
    for k = 1:size(cmbindx,1)
      [tmp, cmbindx(k,:)] = match_str(channelcmb(k,:)', freq.label);
      if ~isempty(intersect(cmbindx(1:(k-1),:), cmbindx(k,:), 'rows'))
        ok(k) = false;
      elseif cmbindx(k,1)==cmbindx(k,2)
        ok(k) = false;
      end
    end
  end
  
  %remove auto-combinations and double occurrences
  cmbindx    = cmbindx(ok,:);
  channelcmb = channelcmb(ok,:);
  
  %do multiple 2x2 factorization efficiently
  if ntim>1,
    for kk = 1:ntim
      [Htmp, Ztmp, Stmp] = sfactorization_wilson2x2(freq.crsspctrm(:,:,:,kk), ...
                               freq.freq, numiteration, tol, cmbindx, fb);
      if kk==1,
        H   = Htmp;
        Z   = Ztmp;
        S   = Stmp;
        H(:,:,2:ntim) = nan;
        Z(:,:,2:ntim) = nan;
        S(:,:,2:ntim) = nan;
      else
        H(:,:,kk) = Htmp;
        Z(:,:,kk) = Ztmp;
        S(:,:,kk) = Stmp;
      end
    end
  else
    % if the number of pairs becomes too big, it seems to slow down quite a
    % bit. try to chunk
    if size(cmbindx,1)>1000
      begchunk = 1:1000:size(cmbindx,1);
      endchunk = [1000:1000:size(cmbindx,1) size(cmbindx,1)];
      H = zeros(4*size(cmbindx,1), numel(freq.freq));
      S = zeros(4*size(cmbindx,1), numel(freq.freq));
      Z = zeros(4*size(cmbindx,1), 1);
      for k = 1:numel(begchunk)
        fprintf('computing factorization of chunck %d/%d\n', k, numel(begchunk));
        [Htmp, Ztmp, Stmp] = sfactorization_wilson2x2(freq.crsspctrm, freq.freq, ...
                                             numiteration, tol, cmbindx(begchunk(k):endchunk(k),:), fb);
                                           
        begix = (k-1)*4000+1;
        endix = min(k*4000, size(cmbindx,1)*4);
        H(begix:endix, :) = Htmp;
        S(begix:endix, :) = Stmp;
        Z(begix:endix, :) = Ztmp;
                                           
      end
    else
      [H, Z, S] = sfactorization_wilson2x2(freq.crsspctrm, freq.freq, ...
                                             numiteration, tol, cmbindx, fb);
    end
  end
  
  %convert crsspctrm accordingly
  siz          = [size(H) 1];
  tmpcrsspctrm = complex(zeros([2 2 siz(1)/4 siz(2:end)]));
  for k = 1:size(cmbindx,1)
    tmpcrsspctrm(:,:,k,:,:) = freq.crsspctrm(cmbindx(k,:),cmbindx(k,:),:,:);   
  end
  freq.crsspctrm = reshape(tmpcrsspctrm, siz);
  
  labelcmb = cell(size(cmbindx,1)*4, 2);
  for k = 1:size(cmbindx,1)
    indx = (k-1)*4 + (1:4);
    labelcmb{indx(1),1} = [channelcmb{k,1},'[',channelcmb{k,1},channelcmb{k,2},']'];
    labelcmb{indx(1),2} = [channelcmb{k,1},'[',channelcmb{k,1},channelcmb{k,2},']'];
    labelcmb{indx(2),1} = [channelcmb{k,2},'[',channelcmb{k,1},channelcmb{k,2},']'];
    labelcmb{indx(2),2} = [channelcmb{k,1},'[',channelcmb{k,1},channelcmb{k,2},']'];
    labelcmb{indx(3),1} = [channelcmb{k,1},'[',channelcmb{k,1},channelcmb{k,2},']'];
    labelcmb{indx(3),2} = [channelcmb{k,2},'[',channelcmb{k,1},channelcmb{k,2},']'];
    labelcmb{indx(4),1} = [channelcmb{k,2},'[',channelcmb{k,1},channelcmb{k,2},']'];
    labelcmb{indx(4),2} = [channelcmb{k,2},'[',channelcmb{k,1},channelcmb{k,2},']'];
  end
elseif ~isempty(block)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % blockwise multivariate stuff
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if ntim>1,
    error('blockwise factorization of tfrs is not yet possible');
  end
  
  %reorder the channel order such that the blocks are ordered
  nblocks = unique(blockindx{2});
  for k = 1:numel(nblocks)
    %b{k} = cfg.blockindx{2}(find(cfg.blockindx{2}==nblocks(k)));
    b{k} = find(blockindx{2}==nblocks(k));
  end
  indx = cat(1,b{:});
  freq.crsspctrm = freq.crsspctrm(indx, indx, :);
  freq.label     = freq.label(indx);
  bindx          = blockindx{2}(indx);

  for k = 1:numel(block)
    sel  = find(ismember(bindx, block{k}));
    Stmp = freq.crsspctrm(sel,sel,:);
    
    % do PCA to avoid zigzags due to numerical issues
    dopca = 1;
    if dopca
      dat     = sum(Stmp,3);
      [u,s,v] = svd(real(dat));
      for m = 1:size(Stmp,3)
        Stmp(:,:,m) = u'*Stmp(:,:,m)*u;
      end
    end
    
    [Htmp, Ztmp, Stmp] = sfactorization_wilson(Stmp, freq.freq, ...
                                                 numiteration, tol, fb);  
    
    % undo PCA
    if dopca
      for m = 1:size(Stmp,3)
        Htmp(:,:,m) = u*Htmp(:,:,m)*u';
        Stmp(:,:,m) = u*Stmp(:,:,m)*u';
      end
      Ztmp = u*Ztmp*u';
    end
                                               
    siz  = [size(Htmp) 1]; siz = [siz(1)*siz(2) siz(3:end)];
    Htmp = reshape(Htmp, siz);
    siz  = [size(Ztmp) 1]; siz = [siz(1)*siz(2) siz(3:end)];
    Ztmp = reshape(Ztmp, siz);
    siz  = [size(Stmp) 1]; siz = [siz(1)*siz(2) siz(3:end)];
    Stmp = reshape(Stmp, siz);
    
    tmpindx = [];
    cmbtmp  = cell(siz(1), 2);
    [tmpindx(:,1), tmpindx(:,2)] = ind2sub(sqrt(siz(1))*[1 1],1:siz(1));
    for kk = 1:size(cmbtmp,1)
      cmbtmp{kk,1} = [freq.label{sel(tmpindx(kk,1))},'[',cat(2,freq.label{sel}),']'];
      cmbtmp{kk,2} = [freq.label{sel(tmpindx(kk,2))},'[',cat(2,freq.label{sel}),']'];
    end
    
    %concatenate
    if k == 1,
      H = Htmp;
      Z = Ztmp;
      S = Stmp;
      labelcmb = cmbtmp;
    else
      H = cat(1,H,Htmp);
      Z = cat(1,Z,Ztmp);
      S = cat(1,S,Stmp);
      labelcmb = cat(1,labelcmb,cmbtmp);
    end 
  end
  
  if strcmp(freq.dimord(1:9), 'chan_chan'),
    freq.dimord = ['chancmb_',freq.dimord(10:end)];
  end
end

% create output
output           = [];
if strcmp(sfmethod, 'multivariate')
  output.dimord    = freq.dimord;
else
  ix = strfind(freq.dimord, 'chan_chan');
  newdimord = [freq.dimord(1:(ix+3)),'cmb',freq.dimord(ix+9:end)];

  %dimtok = tokenize(freq.dimord, '_');
  %chdim  = 0;
  %newdimord = '';
  %for k = 1:numel(dimtok)
  %  if strcmp(dimtok{k}, 'chan'), chdim = chdim+1; end
  %  if chdim==1,
  %    newdimord = [newdimord,'_chancmb'];
  %  elseif chdim
  %    newdimord = [newdimord, '_', dimtok{k}];
  %  end
  %end
  output.dimord    = newdimord;
end
output.label     = freq.label;
output.freq      = freq.freq;
output.crsspctrm = S;
output.transfer  = H;
output.noisecov  = Z;
if isfield(freq, 'time'),      output.time      = freq.time;      end
if isfield(freq, 'cumtapcnt'), output.cumtapcnt = freq.cumtapcnt; end
if isfield(freq, 'cumsucmtn'), output.cumsumcnt = freq.cumsumcnt; end
if exist('labelcmb', 'var') && ~isempty(labelcmb)
  output.labelcmb = labelcmb;
  output          = rmfield(output, 'label');
end
