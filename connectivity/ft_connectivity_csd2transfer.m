function [output] = ft_connectivity_csd2transfer(freq, varargin)

% FT_CONNECTIVITY_CSD2TRANSFER computes the transfer-function from frequency domain
% data using the Wilson-Burg algorithm. The transfer function can be used for the
% computation of directional measures of connectivity, such as Granger causality,
% partial directed coherence, or directed transfer functions.
%
% Use as
%   [output] = ft_connectivity_csd2transfer(freq, ...)
%
% The input variable freq should be a FieldTrip data structure containing frequency
% domain data containing the cross-spectral density computed between all pairs of
% channels, thus containing a 'dimord' of 'chan_chan_freq(_time)'.
%
% Additional optional input arguments come as key-value pairs:
%   numiteration = scalar value (default: 100) the number of iterations
%   channelcmb   = Nx2 cell-array listing the channel pairs for the spectral
%                    factorization. If not defined or empty (default), a
%                    full multivariate factorization is performed, otherwise
%                    a multiple pairwise factorization is done.
%   tol          = scalar value (default: 1e-18) tolerance limit truncating
%                    the iterations
%   sfmethod     = 'multivariate', or 'bivariate'
%   stabilityfix = false, or true. zigzag-reduction by means of tapering of the
%                    intermediate time domain representation when computing the
%                    plusoperator
%
% The code for the Wilson-Burg algorithm has been very generously provided by Dr.
% Mukesh Dhamala, and Prof. Mingzhou Ding and his group, and has been adjusted for
% efficiency. If you use this code for studying directed interactions, please cite
% the following references:
%   - M.Dhamala, R.Rangarajan, M.Ding, Physical Review Letters 100, 018701 (2008).
%   - M.Dhamala, R.Rangarajan, M.Ding, Neuroimage 41, 354 (2008).
%
% See also FT_CONNECTIVITYANALYSIS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Undocumented options:
%   block
%   blockindx
%   svd
%   conditional
%   init
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

sfmethod         = ft_getopt(varargin, 'sfmethod',     'multivariate');
numiteration     = ft_getopt(varargin, 'numiteration', 100);
fb               = ft_getopt(varargin, 'feedback',     'none');
checkconvergence = ft_getopt(varargin, 'checkconvergence', true);
init             = ft_getopt(varargin, 'init',         'chol');
channelcmb       = ft_getopt(varargin, 'channelcmb',   {});
channeltriplet   = ft_getopt(varargin, 'channeltriplet', {});
block            = ft_getopt(varargin, 'block',        []);
tol              = ft_getopt(varargin, 'tol',          1e-18);

dosvd            = istrue(ft_getopt(varargin, 'svd',          'no'));
doconditional    = istrue(ft_getopt(varargin, 'conditional',  0));
stabilityfix     = istrue(ft_getopt(varargin, 'stabilityfix', false));
doblock          = isstruct(block) || doconditional;

if doconditional && isempty(block)
  % create the default block struct-array
  block = struct('name','','label',{});
  for k = 1:numel(freq.label)
    block(k).name  = freq.label{k};
    block(k).label = freq.label(k);
  end
end

% check whether input data is valid, it should contain crsspctrm with label
% (so no linearly indexed labelcmb)
freq = ft_checkdata(freq, 'datatype', 'freq');
if ~isfield(freq, 'crsspctrm') || ~isfield(freq, 'label')
  ft_error('the input data does not contain cross-spectral density data in the supported format');
end

hasrpt = ~isempty(strfind(freq.dimord, 'rpt'));
if hasrpt
  nrpt = numel(freq.cumtapcnt);
else
  nrpt = 1;
end

if ~isempty(channelcmb) && ~(strcmp(channelcmb{1}, 'all') && strcmp(channelcmb{2}, 'all')) && strcmp(sfmethod, 'bivariate')
  % enforce sfmethod to be bivariate
  sfmethod = 'bivariate';
end

if strcmp(sfmethod, 'bivariate')
  % this method allows the {'all' 'all'} channelcmb
  if isempty(channelcmb)
    channelcmb = {'all' 'all'};
  end
elseif strcmp(sfmethod, 'trivariate')
  % this method requires an explicit definition of the triplets (an {'all'
  % 'all' 'all'} scenario might become computationally heavy quite rapidly
  if isempty(channeltriplet)
    ft_error('triplet wise factorization requires an explicit specification of the triplets');
  end
  
  % create an Ntriplet x 3 index matrix, to be used below
  [dum,i2] = match_str(channeltriplet(:,1),freq.label);
  cmbindx(:,1) = i2;
  [dum,i2] = match_str(channeltriplet(:,2),freq.label);
  cmbindx(:,2) = i2;
  [dum,i2] = match_str(channeltriplet(:,3),freq.label);
  cmbindx(:,3) = i2;
  
elseif strcmp(sfmethod, 'bivariate_conditional')
  % this method requires a non-empty channelcmb {Nx2} with a corresponding
  % channelcnd{Nx1}. The elements in channelcnd can be cell-arrays
  if isempty(channelcmb) || size(channelcmb,2) ~= 3
    error('a decomposition that aims at bivariate granger, conditioned on a third channel, requires a channelcmb consisting of 3 columns');
  end
  
  % this method requires a channeltriplet
  channeltriplet = cell(0,3);
  
  % ensure channelcnd to be a cell array of cell arrays
  for k = 1:size(channelcmb,1)
    if ~iscell(channelcmb{k,3})
      channelcmb{k,3} = {channelcmb{k,3}};
    end
  end
  
  for k = 1:size(channelcmb,1)
    for m = 1:numel(channelcmb{k,3})
      channeltriplet(end+1,1:2) = channelcmb(k,1:2);
      channeltriplet{end, 3}    = channelcmb{k,  3}{m};
    end
  end
  
  % remove double occurrences, unique does not work for cell arrays with
  % the argument 'rows'
  tmp     = channeltriplet;
  utmp    = unique(tmp(:));
  ok      = true(size(tmp,1),1);
  tmpindx = zeros(size(tmp));
  for k = 1:size(tmpindx,1)
    [dum, tmpindx(k,:)] = match_str(tmp(k,:)', utmp);
    tmptmpindx = tmpindx(k,:);
    tmptmpindx = tmptmpindx([1 2 3;1 3 2;2 1 3;2 3 1;3 1 2;3 2 1]);
    if ~isempty(intersect(tmpindx(1:(k-1),:), tmptmpindx, 'rows'))
      ok(k) = false;
    elseif tmpindx(k,1)==tmpindx(k,2)
      ok(k) = false;
    end
  end
  channeltriplet = channeltriplet(ok,:);
  %
  %   tmp = cell(0,3);
  %   while ~isempty(channeltriplet)
  %     tmp = cat(1, tmp, channeltriplet(end,:));
  %     channeltriplet(strcmp(channeltriplet(:,1),channeltriplet{end,1})&...
  %                    strcmp(channeltriplet(:,2),channeltriplet{end,2})&...
  %                    strcmp(channeltriplet(:,3),channeltriplet{end,3}),:) = [];
  %   end
  %   channeltriplet = tmp;
  %
end

if ~isempty(channelcmb) && numel(freq.label)>1 && strncmp(sfmethod, 'bivariate', 9)
  if numel(channelcmb)==2 && strcmp(channelcmb{1},'all') && strcmp(channelcmb{2}, 'all')
    % create an Npair x 2 indexing matrix, to be used below
    [cmbindx(:,1), cmbindx(:,2)] = find(tril(ones(numel(freq.label)),-1));
    ok = true(size(cmbindx,1),1);
  end
  channelcmb = ft_channelcombination(channelcmb(:,1:2), freq.label);
end

if ~isempty(block)
  % sanity check 1
  if ~isstruct(block)
    ft_error('block should be a struct-array');
  end
  % sanity check 2
  if strcmp(sfmethod, 'bivariate')
    ft_error('when block is specified, it is not OK to do bivariate decomposition');
  end
end

if isfield(freq, 'time')
  ntim = numel(freq.time);
else
  ntim = 1;
end

siz = size(freq.crsspctrm);
if ntim==1
  siz = [siz 1]; % add dummy dimensionality for time axis
end

if strcmp(sfmethod, 'bivariate')
  list1 = char(channelcmb(:,1)); siz1 = size(list1,2);
  list2 = char(channelcmb(:,2)); siz2 = size(list2,2);
  ix    = 1:min(siz1,siz2);
  npair = sum(sum(list1(:,ix)~=list2(:,ix),2)~=0);
  fprintf('computing pairwise non-parametric spectral factorization on %d channel pairs\n', npair);%size(channelcmb,1) - numel(unique(channelcmb(:))));
elseif strcmp(sfmethod, 'bivariate_conditional')
  % no text
elseif strcmp(sfmethod, 'multivariate')
  fprintf('computing multivariate non-parametric spectral factorization on %d channels\n', numel(freq.label));
elseif strcmp(sfmethod, 'trivariate')
  fprintf('computing tripletwise non-parametric spectral factorization\n');
else
  ft_error('unknown sfmethod %s', sfmethod);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computations start here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(sfmethod, 'multivariate') && nrpt==1 && ~doconditional
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
    
    if any(isnan(tmp(:)))
      Htmp = nan;
      Ztmp = nan;
      Stmp = nan;
    else
      [Htmp, Ztmp, Stmp] = sfactorization_wilson(tmp, freq.freq, ...
        numiteration, tol, fb, init, checkconvergence, stabilityfix);
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
  
elseif strcmp(sfmethod, 'multivariate') && nrpt==1 && doblock
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % blockwise multivariate stuff
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if ntim>1
    ft_error('blockwise factorization of tfrs is not yet possible');
  end
  
  % create a blockindx array that assigns each channel to a block
  blockindx = zeros(numel(freq.label),1);
  for k = 1:numel(block)
    % FIXME make robust against double occurrences
    [ix, iy]      = match_str(freq.label, block(k).label);
    blockindx(ix) = k;
  end
  
  if doconditional
    % multiple factorizations are needed:
    % -all blocks together
    % -all pairs of blocks (
    % -all leave-one-out blocks
    blocks  = unique(blockindx);
    nblocks = numel(blocks);
    
    tmp1 = nchoosek(blocks,2);
    if nblocks>3
      tmp2 = nchoosek(blocks,nblocks-1);
    else
      tmp2 = [];
    end
    
    nfact = size(tmp1,1)+size(tmp2,1)+1;
    maxnfact = 500;
    if nfact>maxnfact
      ft_error('at present the number of factorizations for conditional granger is set to 500');
    end
    
    factorizations = cell(nfact,1);
    factorizations{1} = blocks(:)';
    for k = 1:size(tmp1,1)
      factorizations{k+1} = tmp1(k,:);
    end
    for k = 1:size(tmp2,1)
      factorizations{k+1+size(tmp1,1)} = tmp2(k,:);
    end
  else
    % only a single factorization is needed
    % this is actually handled above
  end
  
  %reorder the channel order such that the blocks are ordered
  for k = 1:nblocks
    %b{k} = cfg.blockindx{2}(find(cfg.blockindx{2}==nblocks(k)));
    b{k} = find(blockindx==blocks(k));
  end
  indx = cat(1,b{:});
  freq.crsspctrm = freq.crsspctrm(indx, indx, :);
  freq.label     = freq.label(indx);
  freq.blockindx = blockindx(indx);
  freq.block     = block;
  
  for k = 1:nfact
    sel  = find(ismember(freq.blockindx, factorizations{k}(:)));
    Stmp = freq.crsspctrm(sel,sel,:);
    
    % do PCA to avoid zigzags due to numerical issues
    dopca = 0;
    if dopca
      dat     = sum(Stmp,3);
      [u,s,v] = svd(real(dat));
      for m = 1:size(Stmp,3)
        Stmp(:,:,m) = u'*Stmp(:,:,m)*u;
      end
    end
    
    [Htmp, Ztmp, Stmp] = sfactorization_wilson(Stmp, freq.freq, ...
      numiteration, tol, fb, init, checkconvergence, stabilityfix);
    
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
      cmbtmp{kk,1} = [freq.label{sel(tmpindx(kk,1))},'[',strjoin(freq.label(sel),','),']'];
      cmbtmp{kk,2} = [freq.label{sel(tmpindx(kk,2))},'[',strjoin(freq.label(sel),','),']'];
    end
    
    %concatenate
    if k == 1
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
  
  if strcmp(freq.dimord(1:9), 'chan_chan')
    freq.dimord = ['chancmb_',freq.dimord(strfind(freq.dimord,'freq'):end)];
  end
  
elseif strcmp(sfmethod, 'multivariate') && nrpt>1 && ~doblock
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
        numiteration, tol, fb, init, checkconvergence, stabilityfix);
      H(k,:,:,:,m) = Htmp;
      Z(k,:,:,m)   = Ztmp;
      S(k,:,:,:,m) = Stmp;
    end
  end
  
elseif strcmp(sfmethod, 'multivariate') && nrpt>1 && doblock && ~doconditional
  % error
  ft_error('single trial estimates and blockwise factorisation is not yet implemented');
  
elseif strcmp(sfmethod, 'bivariate')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % pairwise factorization resulting in linearly indexed transfer functions
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %convert list of channel labels into indices
  if ~exist('cmbindx', 'var')
    cmbindx     = zeros(size(channelcmb));
    ok          = true(size(cmbindx,1), 1);
    for k = 1:size(cmbindx,1)
      [dum, cmbindx(k,:)] = match_str(channelcmb(k,:)', freq.label);
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
  if ntim>1
    for kk = 1:ntim
      [Htmp, Ztmp, Stmp] = sfactorization_wilson2x2(freq.crsspctrm(:,:,:,kk), ...
        freq.freq, numiteration, tol, cmbindx, fb, init, checkconvergence, stabilityfix);
      if kk==1
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
    nperchunk = 2000;
    if size(cmbindx,1)>nperchunk
      begchunk = 1:nperchunk:size(cmbindx,1);
      endchunk = [nperchunk:nperchunk:size(cmbindx,1) size(cmbindx,1)];
      H = zeros(4*size(cmbindx,1), numel(freq.freq));
      S = zeros(4*size(cmbindx,1), numel(freq.freq));
      Z = zeros(4*size(cmbindx,1), 1);
      for k = 1:numel(begchunk)
        fprintf('computing factorization of chunck %d/%d\n', k, numel(begchunk));
        [Htmp, Ztmp, Stmp] = sfactorization_wilson2x2(freq.crsspctrm, freq.freq, ...
          numiteration, tol, cmbindx(begchunk(k):endchunk(k),:), fb, init, checkconvergence, stabilityfix);
        
        begix = (k-1)*nperchunk*4+1;
        endix = min(k*nperchunk*4, size(cmbindx,1)*4);
        H(begix:endix, :) = Htmp;
        S(begix:endix, :) = Stmp;
        Z(begix:endix, :) = Ztmp;
        
      end
    else
      [H, Z, S] = sfactorization_wilson2x2(freq.crsspctrm, freq.freq, ...
        numiteration, tol, cmbindx, fb, init, checkconvergence, stabilityfix);
    end
  end
  
  labelcmb = cell(size(cmbindx,1)*4, 2);
  for k = 1:size(cmbindx,1)
    duplet = sprintf('%s,%s',channelcmb{k,1},channelcmb{k,2});
    indx   = (k-1)*4 + (1:4);
    labelcmb{indx(1),1} = [channelcmb{k,1},'[',duplet,']'];
    labelcmb{indx(1),2} = [channelcmb{k,1},'[',duplet,']'];
    labelcmb{indx(2),1} = [channelcmb{k,2},'[',duplet,']'];
    labelcmb{indx(2),2} = [channelcmb{k,1},'[',duplet,']'];
    labelcmb{indx(3),1} = [channelcmb{k,1},'[',duplet,']'];
    labelcmb{indx(3),2} = [channelcmb{k,2},'[',duplet,']'];
    labelcmb{indx(4),1} = [channelcmb{k,2},'[',duplet,']'];
    labelcmb{indx(4),2} = [channelcmb{k,2},'[',duplet,']'];
  end
  
elseif strcmp(sfmethod, 'bivariate_conditional')
  % recursively call this function, and concatenate the output
  
  optarg = varargin;
  sel    = strcmp(optarg, 'bivariate_conditional');
  optarg{sel} = 'trivariate';
  sel    = find(strcmp(optarg, 'channeltriplet'));
  if ~isempty(sel)
    optarg{sel+1} = channeltriplet;
  else
    optarg = cat(2, optarg, {'channeltriplet', channeltriplet});
  end
  out3   = ft_connectivity_csd2transfer(freq, optarg{:});
  
  % get the required duplets
  optarg = varargin;
  sel    = strcmp(optarg, 'bivariate_conditional');
  optarg{sel} = 'bivariate';
  sel    = find(strcmp(optarg, 'channelcmb'));
  
  % ensure all the necessary duplets, for bidirectionality, but exclude
  % double occurrences
  tmp     = [optarg{sel+1}(:,[1 3]);optarg{sel+1}(:,[2 3])];
  utmp    = unique(tmp(:));
  ok      = true(size(tmp,1),1);
  tmpindx = zeros(size(tmp));
  for k = 1:size(tmpindx,1)
    [dum, tmpindx(k,:)] = match_str(tmp(k,:)', utmp);
    if ~isempty(intersect(tmpindx(1:(k-1),:), [tmpindx(k,:);tmpindx(k,[2 1])], 'rows'))
      ok(k) = false;
    elseif tmpindx(k,1)==tmpindx(k,2)
      ok(k) = false;
    end
  end
  optarg{sel+1} = tmp(ok,:);
  
  out2   = ft_connectivity_csd2transfer(freq, optarg{:});
  
  H = cat(1, out3.transfer,  out2.transfer);
  Z = cat(1, out3.noisecov,  out2.noisecov);
  S = cat(1, out3.crsspctrm, out2.crsspctrm);
  labelcmb = cat(1, out3.labelcmb, out2.labelcmb);
  
elseif strcmp(sfmethod, 'trivariate')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % tripletwise factorization resulting in linearly indexed transfer functions
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %FIXME remove auto-combinations and double occurrences
  
  %do multiple 3x3 factorization efficiently
  if ntim>1
    for kk = 1:ntim
      [Htmp, Ztmp, Stmp] = sfactorization_wilson3x3(freq.crsspctrm(:,:,:,kk), ...
        freq.freq, numiteration, tol, cmbindx, fb, init, checkconvergence, stabilityfix);
      if kk==1
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
    nperchunk = 2000;
    if size(cmbindx,1)>nperchunk
      begchunk = 1:nperchunk:size(cmbindx,1);
      endchunk = [nperchunk:nperchunk:size(cmbindx,1) size(cmbindx,1)];
      H = zeros(9*size(cmbindx,1), numel(freq.freq));
      S = zeros(9*size(cmbindx,1), numel(freq.freq));
      Z = zeros(9*size(cmbindx,1), 1);
      for k = 1:numel(begchunk)
        fprintf('computing factorization of chunck %d/%d\n', k, numel(begchunk));
        [Htmp, Ztmp, Stmp] = sfactorization_wilson3x3(freq.crsspctrm, freq.freq, ...
          numiteration, tol, cmbindx(begchunk(k):endchunk(k),:), fb, init, checkconvergence, stabilityfix);
        
        begix = (k-1)*nperchunk*9+1;
        endix = min(k*nperchunk*9, size(cmbindx,1)*9);
        H(begix:endix, :) = Htmp;
        S(begix:endix, :) = Stmp;
        Z(begix:endix, :) = Ztmp;
        
      end
    else
      [H, Z, S] = sfactorization_wilson3x3(freq.crsspctrm, freq.freq, ...
        numiteration, tol, cmbindx, fb, init, checkconvergence, stabilityfix);
    end
  end
  
  labelcmb = cell(size(cmbindx,1)*9, 2);
  for k = 1:size(cmbindx,1)
    triplet = sprintf('%s,%s,%s',channeltriplet{k,1},channeltriplet{k,2},channeltriplet{k,3});
    indx    = (k-1)*9 + (1:9);
    
    labelcmb{indx(1),1} = [channeltriplet{k,1},'[',triplet,']'];
    labelcmb{indx(1),2} = [channeltriplet{k,1},'[',triplet,']'];
    labelcmb{indx(2),1} = [channeltriplet{k,2},'[',triplet,']'];
    labelcmb{indx(2),2} = [channeltriplet{k,1},'[',triplet,']'];
    labelcmb{indx(3),1} = [channeltriplet{k,3},'[',triplet,']'];
    labelcmb{indx(3),2} = [channeltriplet{k,1},'[',triplet,']'];
    labelcmb{indx(4),1} = [channeltriplet{k,1},'[',triplet,']'];
    labelcmb{indx(4),2} = [channeltriplet{k,2},'[',triplet,']'];
    labelcmb{indx(5),1} = [channeltriplet{k,2},'[',triplet,']'];
    labelcmb{indx(5),2} = [channeltriplet{k,2},'[',triplet,']'];
    labelcmb{indx(6),1} = [channeltriplet{k,3},'[',triplet,']'];
    labelcmb{indx(6),2} = [channeltriplet{k,2},'[',triplet,']'];
    labelcmb{indx(7),1} = [channeltriplet{k,1},'[',triplet,']'];
    labelcmb{indx(7),2} = [channeltriplet{k,3},'[',triplet,']'];
    labelcmb{indx(8),1} = [channeltriplet{k,2},'[',triplet,']'];
    labelcmb{indx(8),2} = [channeltriplet{k,3},'[',triplet,']'];
    labelcmb{indx(9),1} = [channeltriplet{k,3},'[',triplet,']'];
    labelcmb{indx(9),2} = [channeltriplet{k,3},'[',triplet,']'];
  end
elseif strcmp(sfmethod, 'bivariate') && nrpt>1
  % error
  ft_error('single trial estimates and linear combination indexing is not implemented');
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% create output structure
%%%%%%%%%%%%%%%%%%%%%%%%%
output = keepfields(freq, {'freq' 'time' 'cumtapcnt' 'cumsumcnt' 'block' 'blockindx'});
output.crsspctrm = S;
output.transfer  = H;
output.noisecov  = Z;
if strcmp(sfmethod, 'multivariate')
  output.dimord    = freq.dimord;
else
  if strncmp(freq.dimord, 'chan_chan', 9)
    freq.dimord = ['chancmb_',freq.dimord(strfind(freq.dimord,'freq'):end)];
  end
  output.dimord    = freq.dimord;
end
if exist('labelcmb', 'var') && ~isempty(labelcmb)
  output.labelcmb = labelcmb;
else
  output.label    = freq.label;
end
