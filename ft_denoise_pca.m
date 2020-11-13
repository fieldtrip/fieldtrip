function data = ft_denoise_pca(cfg, varargin)

% FT_DENOISE_PCA performs a principal component analysis (PCA) on specified reference
% channels and subtracts the projection of the data of interest onto this orthogonal
% basis from the data of interest. This is the algorithm which is applied by 4D to
% compute noise cancellation weights on a dataset of interest. This function has been
% designed for 4D MEG data, but can also be applied to data from other MEG systems.
%
% Use as
%   [dataout] = ft_denoise_pca(cfg, data)
% or as
%   [dataout] = ft_denoise_pca(cfg, data, refdata)
% where "data" is a raw data structure that was obtained with FT_PREPROCESSING. If
% you specify the additional input "refdata", the specified reference channels for
% the regression will be taken from this second data structure. This can be useful
% when reference-channel specific preprocessing needs to be done (e.g. low-pass
% filtering).
%
% The output structure dataout contains the denoised data in a format that is
% consistent with the output of FT_PREPROCESSING.
%
% The configuration should contain
%   cfg.refchannel = the channels used as reference signal (default = 'MEGREF')
%   cfg.channel    = the channels to be denoised (default = 'MEG')
%   cfg.truncate   = optional truncation of the singular value spectrum (default = 'no')
%   cfg.zscore     = standardise reference data prior to PCA (default = 'no')
%   cfg.pertrial   = 'no' (default) or 'yes'. Regress out the references on a per trial basis
%   cfg.trials     = list of trials that are used (default = 'all')
%   cfg.updatesens = 'no' or 'yes' (default = 'yes')
%
% if cfg.truncate is integer n > 1, n will be the number of singular values kept.
% if 0 < cfg.truncate < 1, the singular value spectrum will be thresholded at the
% fraction cfg.truncate of the largest singular value.
%
% See also FT_PREPROCESSING, FT_DENOISE_SYNTHETIC

% Undocumented cfg-option: cfg.pca the output structure of an earlier call
% to the function. Can be used regress out the reference channels from
% another data set.

% Copyright (c) 2008-2009, Jan-Mathijs Schoffelen, CCNi Glasgow
% Copyright (c) 2010-2011, Jan-Mathijs Schoffelen, DCCN Nijmegen
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble provenance varargin
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'raw');
end

% set the defaults
cfg.refchannel = ft_getopt(cfg, 'refchannel', 'MEGREF');
cfg.channel    = ft_getopt(cfg, 'channel',    'MEG');
cfg.truncate   = ft_getopt(cfg, 'truncate',   'no');
cfg.zscore     = ft_getopt(cfg, 'zscore',     'no');
cfg.trials     = ft_getopt(cfg, 'trials',     'all', 1);
cfg.pertrial   = ft_getopt(cfg, 'pertrial',   'no');
cfg.feedback   = ft_getopt(cfg, 'feedback',   'none');
cfg.updatesens = ft_getopt(cfg, 'updatesens', 'yes');


if istrue(cfg.pertrial)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % iterate over trials
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  tmpcfg  = keepfields(cfg, {'trials', 'showcallinfo'});
  % select trials of interest
  for i=1:numel(varargin)
    varargin{i}        = ft_selectdata(tmpcfg, varargin{i});
    [dum, varargin{i}] = rollback_provenance(cfg, varargin{i});
    if i==1
      cfg = dum;
    end
  end

  tmp             = cell(numel(varargin{1}.trial),1);
  tmpcfg          = cfg;
  tmpcfg.pertrial = 'no';
  for k = 1:numel(varargin{1}.trial)
    tmpcfg.trials   = k;    % select a single trial
    tmp{k}          = ft_denoise_pca(tmpcfg, varargin{:});
    [dum, tmp{k}]   = rollback_provenance(tmpcfg, tmp{k});
  end
  data = ft_appenddata([], tmp{:});
  [cfg, data] = rollback_provenance(cfg, data);

else
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compute it for the data concatenated over all trials
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  computeweights = ~isfield(cfg, 'pca');

  % select trials of interest
  tmpcfg  = keepfields(cfg, {'trials', 'showcallinfo'});
  if length(varargin)==1
    % channel data and reference channel data are in 1 data structure
    megchan = ft_channelselection(cfg.channel,    varargin{1}.label);
    refchan = ft_channelselection(cfg.refchannel, varargin{1}.label);
      
    tmpcfg.channel = refchan;
    refdata        = ft_selectdata(tmpcfg, varargin{1});
    [dum,refdata]  = rollback_provenance(cfg, refdata);
    tmpcfg.channel = megchan;
    data           = ft_selectdata(tmpcfg, varargin{1});
    [cfg, data]    = rollback_provenance(cfg, data);
  
  else
    % channel data and reference channel data are in 2 data structures
    megchan = ft_channelselection(cfg.channel,    varargin{1}.label);
    refchan = ft_channelselection(cfg.refchannel, varargin{2}.label);
    
    % throw a warning if some of the specified reference channels are also
    % in the first data argument
    if ~isempty(ft_channelselection(cfg.refchannel, varargin{1}.label))
      ft_warning('some of the specified reference channels are also present in the first data argument, this information will not be used for the cleaning of the data');
    end
    
    tmpcfg.channel = refchan;
    refdata        = ft_selectdata(tmpcfg, varargin{2});
    [dum, refdata] = rollback_provenance(cfg, refdata);
    tmpcfg.channel = megchan;
    data           = ft_selectdata(tmpcfg, varargin{1});
    [cfg, data]    = rollback_provenance(cfg, data);
    
  end
  
  refchan = ft_channelselection(cfg.refchannel, refdata.label);
  refindx = match_str(refdata.label, refchan);
  megchan = ft_channelselection(cfg.channel, data.label);
  megindx = match_str(data.label, megchan);

  nref = length(refindx);
  ntrl = length(data.trial);

  if ischar(cfg.truncate) && strcmp(cfg.truncate, 'no')
    cfg.truncate = length(refindx);
  elseif ischar(cfg.truncate) || (cfg.truncate>1 && cfg.truncate/round(cfg.truncate)~=1) || cfg.truncate>length(refindx)
    ft_error('cfg.truncate should be either ''no'', an integer number <= the number of references, or a number between 0 and 1');
    % FIXME the default truncation applied by 4D is 1x10^-8
  end

  % compute and remove mean from data
  fprintf('removing the mean from the channel data and reference channel data\n');
  m             = cellmean(data.trial,       2);
  data.trial    = cellvecadd(data.trial,    -m);
  m             = cellmean(refdata.trial,    2);
  refdata.trial = cellvecadd(refdata.trial, -m);

  % compute std of data before the regression
  stdpre = cellstd(data.trial, 2);

  if computeweights

    % zscore
    if strcmp(cfg.zscore, 'yes')
      fprintf('zscoring the reference channel data\n');
      [refdata.trial, sdref] = cellzscore(refdata.trial, 2, 0); %forced demeaned already
    else
      sdref = ones(nref, 1);
    end

    % compute covariance of refchannels and do svd
    fprintf('performing pca on the reference channel data\n');
    crefdat = cellcov(refdata.trial, [], 2, 0);
    [u,s,v] = svd(crefdat);

    % determine the truncation and rotation
    if cfg.truncate<1
      % keep all singular vectors with singular values >= cfg.truncate*s(1,1)
      s1   = s./max(s(:));
      keep = find(diag(s1)>cfg.truncate);
    else
      keep = 1:cfg.truncate;
    end
    fprintf('keeping %d out of %d components\n',numel(keep),size(u,2));
    rotmat = u(:, keep)';

    % rotate the refdata
    fprintf('projecting the reference data onto the pca-subspace\n');
    refdata.trial = cellfun(@mtimes, repmat({rotmat}, 1, ntrl), refdata.trial, 'UniformOutput', 0);

    % project megdata onto the orthogonal basis
    fprintf('computing the regression weights\n');
    nom   = cellcov(data.trial,    refdata.trial, 2, 0);
    denom = cellcov(refdata.trial, [],            2, 0);
    rw    = (pinv(denom)*nom')';

    % subtract projected data
    fprintf('subtracting the reference channel data from the channel data\n');
    for k = 1:ntrl
      data.trial{k} = data.trial{k} - rw*refdata.trial{k};
    end

    % rotate back and 'unscale'
    pca.w        = rw*rotmat*diag(1./sdref);
    pca.label    = data.label;
    pca.reflabel = refdata.label;
    pca.rotmat   = rotmat;
    cfg.pca      = pca;

  else
    fprintf('applying precomputed weights to the data\n');
    % check whether the weight table contains the specified references
    % ensure the ordering of the meg-data to be consistent with the weights
    % ensure the ordering of the ref-data to be consistent with the weights

    [i1,i2] = match_str(refchan, cfg.pca.reflabel);
    [i3,i4] = match_str(megchan, cfg.pca.label);
    if length(i2)~=length(cfg.pca.reflabel)
      ft_error('you specified fewer references to use as there are in the precomputed weight table');
    end

    refindx = refindx(i1);
    megindx = megindx(i3);
    cfg.pca.w = cfg.pca.w(i4,i2);
    cfg.pca.label   = cfg.pca.label(i4);
    cfg.pca.reflabel= cfg.pca.reflabel(i2);
    if isfield(cfg.pca, 'rotmat')
      cfg.pca = rmfield(cfg.pca, 'rotmat'); % dont know
    end

    for k = 1:ntrl
      data.trial{k} = data.trial{k} - cfg.pca.w*refdata.trial{k};
    end
    pca = cfg.pca;

  end

  % compute std of data after
  stdpst = cellstd(data.trial, 2);

  % demean FIXME is this needed
  m          = cellmean(data.trial, 2);
  data.trial = cellvecadd(data.trial, -m);

  if isfield(data, 'grad')
    sensfield = 'grad';
  elseif isfield(data, 'elec')
    sensfield = 'elec';
  elseif isfield(data, 'opto')
    sensfield = 'opto';
  else
    sensfield = [];
  end

  % apply the linear projection also to the sensor description
  if ~isempty(sensfield)
    if  strcmp(cfg.updatesens, 'yes')
      fprintf('also applying the weights to the %s structure\n', sensfield);

      montage     = [];
      labelnew    = pca.label;

      % add columns of refchannels not yet present in labelnew
      % [id, i1]  = setdiff(pca.reflabel, labelnew);
      % labelold  = [labelnew; pca.reflabel(sort(i1))];
      labelold  = data.grad.label;
      nlabelold = length(labelold);

      % start with identity
      montage.tra = eye(nlabelold);

      % subtract weights
      [i1, i2]  = match_str(labelold, pca.reflabel);
      [i3, i4]  = match_str(labelold, pca.label);
      montage.tra(i3,i1) = montage.tra(i3,i1) - pca.w(i4,i2);
      montage.labelold  = labelold;
      montage.labelnew  = labelold;

      data.(sensfield) = ft_apply_montage(data.(sensfield), montage, 'keepunused', 'yes', 'balancename', 'pca');

      % order the fields
      fnames = fieldnames(data.(sensfield).balance);
      tmp    = false(1,numel(fnames));
      for k = 1:numel(fnames)
        tmp(k) = isstruct(data.(sensfield).balance.(fnames{k}));
      end
      [tmp, ix] = sort(tmp, 'descend');
      data.grad.balance = orderfields(data.(sensfield).balance, fnames(ix));

    else
      fprintf('not applying the weights to the %s structure\n', sensfield);
    end
  end % if sensfield

end % if pertrial

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   varargin
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data

%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%

%-----cellcov
function [c] = cellcov(x, y, dim, flag)

% [C] = CELLCOV(X, DIM) computes the covariance, across all cells in x along
% the dimension dim. When there are three inputs, covariance is computed between
% all cells in x and y
%
% X (and Y) should be linear cell-array(s) of matrices for which the size in at
% least one of the dimensions should be the same for all cells

if nargin==2
  flag = 1;
  dim  = y;
  y    = [];
elseif nargin==3
  flag = 1;
end

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1)
  ft_error('incorrect input for cellmean');
end

if nargin==1
  scx1 = cellfun('size', x, 1);
  scx2 = cellfun('size', x, 2);
  if     all(scx2==scx2(1)), dim = 2; %let second dimension prevail
  elseif all(scx1==scx1(1)), dim = 1;
  else   ft_error('no dimension to compute covariance for');
  end
end

if flag
  mx   = cellmean(x, 2);
  x    = cellvecadd(x, -mx);
  if ~isempty(y)
    my = cellmean(y, 2);
    y  = cellvecadd(y, -my);
  end
end

nx   = max(nx);
nsmp = cellfun('size', x, dim);
if isempty(y)
  csmp = cellfun(@covc, x, repmat({dim},1,nx), 'UniformOutput', 0);
else
  csmp = cellfun(@covc, x, y, repmat({dim},1,nx), 'UniformOutput', 0);
end
nc   = size(csmp{1});
c    = sum(reshape(cell2mat(csmp), [nc(1) nc(2) nx]), 3)./sum(nsmp);

function [c] = covc(x, y, dim)

if nargin==2
  dim = y;
  y   = x;
end

if dim==1
  c = x'*y;
elseif dim==2
  c = x*y';
end

%-----cellmean
function [m] = cellmean(x, dim)

% [M] = CELLMEAN(X, DIM) computes the mean, across all cells in x along
% the dimension dim.
%
% X should be an linear cell-array of matrices for which the size in at
% least one of the dimensions should be the same for all cells

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1)
  ft_error('incorrect input for cellmean');
end

if nargin==1
  scx1 = cellfun('size', x, 1);
  scx2 = cellfun('size', x, 2);
  if     all(scx2==scx2(1)), dim = 2; %let second dimension prevail
  elseif all(scx1==scx1(1)), dim = 1;
  else   ft_error('no dimension to compute mean for');
  end
end

nx   = max(nx);
nsmp = cellfun('size', x, dim);
ssmp = cellfun(@sum,   x, repmat({dim},1,nx), 'UniformOutput', 0);
m    = sum(cell2mat(ssmp), dim)./sum(nsmp);

%-----cellstd
function [sd] = cellstd(x, dim, flag)

% [M] = CELLSTD(X, DIM, FLAG) computes the standard deviation, across all cells in x along
% the dimension dim, normalising by the total number of samples
%
% X should be an linear cell-array of matrices for which the size in at
% least one of the dimensions should be the same for all cells. If flag==1, the mean will
% be subtracted first (default behavior, but to save time on already demeaned data, it
% can be set to 0).

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1)
  ft_error('incorrect input for cellstd');
end

if nargin<2
  scx1 = cellfun('size', x, 1);
  scx2 = cellfun('size', x, 2);
  if     all(scx2==scx2(1)), dim = 2; %let second dimension prevail
  elseif all(scx1==scx1(1)), dim = 1;
  else   ft_error('no dimension to compute mean for');
  end
elseif nargin==2
  flag = 1;
end

if flag
  m    = cellmean(x, dim);
  x    = cellvecadd(x, -m);
end

nx   = max(nx);
nsmp = cellfun('size', x, dim);
ssmp = cellfun(@sumsq,   x, repmat({dim},1,nx), 'UniformOutput', 0);
sd   = sqrt(sum(cell2mat(ssmp), dim)./sum(nsmp));

function [s] = sumsq(x, dim)

s = sum(x.^2, dim);

%-----cellvecadd
function [y] = cellvecadd(x, v)

% [Y]= CELLVECADD(X, V) - add vector to all rows or columns of each matrix
% in cell-array X

% check once and for all to save time
persistent bsxfun_exists;
if isempty(bsxfun_exists)
  bsxfun_exists=exist('bsxfun', 'builtin');
  if ~bsxfun_exists
    ft_error('bsxfun not found.');
  end
end

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1)
  ft_error('incorrect input for cellmean');
end

if ~iscell(v)
  v = repmat({v}, nx);
end

y  = cellfun(@bsxfun, repmat({@plus}, nx), x, v, 'UniformOutput', 0);

%-----cellvecmult
function [y] = cellvecmult(x, v)

% [Y]= CELLVECMULT(X, V) - multiply vectors in cell-array V
% to all rows or columns of each matrix in cell-array X
% V can be a vector or a cell-array of vectors

% check once and for all to save time
persistent bsxfun_exists;
if isempty(bsxfun_exists)
  bsxfun_exists=exist('bsxfun', 'builtin');
  if ~bsxfun_exists
    ft_error('bsxfun not found.');
  end
end

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1)
  ft_error('incorrect input for cellmean');
end

if ~iscell(v)
  v = repmat({v}, nx);
end

sx1 = cellfun('size', x, 1);
sx2 = cellfun('size', x, 2);
sv1 = cellfun('size', v, 1);
sv2 = cellfun('size', v, 2);
if     all(sx1==sv1) && all(sv2==1)
elseif all(sx2==sv2) && all(sv1==1)
elseif all(sv1==1)   && all(sv2==1)
else   ft_error('inconsistent input');
end

y  = cellfun(@bsxfun, repmat({@times}, nx), x, v, 'UniformOutput', 0);

%-----cellzscore
function [z, sd, m] = cellzscore(x, dim, flag)

% [Z, SD] = CELLZSCORE(X, DIM, FLAG) computes the zscore, across all cells in x along
% the dimension dim, normalising by the total number of samples
%
% X should be an linear cell-array of matrices for which the size in at
% least one of the dimensions should be the same for all cells. If flag==1, the mean will
% be subtracted first (default behavior, but to save time on already demeaned data, it
% can be set to 0). SD is a vector containing the standard deviations, used for the normalisation.

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1)
  ft_error('incorrect input for cellstd');
end

if nargin<2
  scx1 = cellfun('size', x, 1);
  scx2 = cellfun('size', x, 2);
  if     all(scx2==scx2(1)), dim = 2; % let second dimension prevail
  elseif all(scx1==scx1(1)), dim = 1;
  else   ft_error('no dimension to compute mean for');
  end
elseif nargin==2
  flag = 1;
end

if flag
  m    = cellmean(x, dim);
  x    = cellvecadd(x, -m);
end

sd   = cellstd(x, dim, 0);
z    = cellvecmult(x, 1./sd);
