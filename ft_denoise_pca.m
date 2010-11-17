function [data, pca, stdpre, stdpst] = ft_denoise_pca(cfg, varargin)

% DENOISE_PCA performs a pca on specified reference channels and subtracts
% the projection of the data of interest onto this orthogonal basis from 
% the data of interest. This is the algorithm which is applied by 4D to 
% compute noise cancellation weights on a dataset of interest.
%
% Use as
%   [data] = ft_denoise_pca(cfg, data) or [data] = ft_denoise_pca(cfg, data, refdata)
%
% The configuration should be according to
%   cfg.refchannel = the channels used as reference signal (default = 'MEGREF')
%   cfg.channel    = the channels to be denoised (default = 'MEG')
%   cfg.truncate   = optional truncation of the singular value spectrum
%   cfg.zscore     = standardise reference data prior to PCA (default = 'no')
%     if cfg.truncate is integer n > 1, n will be the number of singular values kept.
%     if 0 < cfg.truncate < 1, the singular value spectrum will be thresholded at the 
%     fraction cfg.truncate of the largest singular value.
%     (default = 'no');

% Copyright (c) 2008-2009, Jan-Mathijs Schoffelen, CCNi Glasgow
% Copyright (c) 2010, Jan-Mathijs Schoffelen, DCCN Nijmegen
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
% $Id: ft_denoise_pca.m 2097 2010-11-10 09:20:18Z roboos $

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

if ~isfield(cfg, 'truncate'),    cfg.truncate   = 'no';     end;
if ~isfield(cfg, 'channel'),     cfg.channel    = 'MEG';    end;
if ~isfield(cfg, 'refchannel'),  cfg.refchannel = 'MEGREF'; end;
if ~isfield(cfg, 'trials'),      cfg.trials     = 'all';    end;
if ~isfield(cfg, 'zscore'),      cfg.zscore     = 'no';     end;
if ~isfield(cfg, 'blc'),         cfg.blc        = 'yes';    end;
if ~isfield(cfg, 'pertrial'),    cfg.pertrial   = 'no';     end

if strcmp(cfg.pertrial, 'yes'),
  tmpcfg          = cfg;
  tmpcfg.pertrial = 'no';
  tmp             = cell(numel(varargin{1}.trial));
  for k = 1:numel(varargin{1}.trial)
    tmpcfg.trials = k;
    tmp{k,1}      = ft_denoise_pca(tmpcfg, varargin{:});
  end
  data = ft_appenddata([], tmp{:});
  return;  
end

computeweights = ~isfield(cfg, 'pca');

%overrule insensible user specifications
%if ~strcmp(cfg.blc, 'yes'),
%  cfg.blc = 'yes';
%end

if length(varargin)==1,
  data    = varargin{1};
  megchan = ft_channelselection(cfg.channel, data.label);
  refchan = ft_channelselection(cfg.refchannel, data.label);
   
  %split data into data and refdata
  tmpcfg  = [];
  tmpcfg.channel = refchan;
  refdata = ft_preprocessing(tmpcfg, data);
  tmpcfg.channel = megchan;
  data    = ft_preprocessing(tmpcfg, data);
else
  data    = varargin{1};
  refdata = varargin{2};
  megchan = ft_channelselection(cfg.channel, data.label);
  refchan = ft_channelselection(cfg.refchannel, refdata.label);
   
  %split data into data and refdata
  tmpcfg  = [];
  tmpcfg.channel = refchan;
  refdata = ft_preprocessing(tmpcfg, refdata);
  tmpcfg.channel = megchan;
  data    = ft_preprocessing(tmpcfg, data);
  
  %FIXME do compatibility check on data vs refdata with respect to dimensions (time-trials)
end

% select trials of interest
if ~isfield(cfg, 'trials'),   cfg.trials = 'all';  end % set the default
if ~strcmp(cfg.trials, 'all')
  fprintf('selecting %d trials\n', length(cfg.trials));
  data    = ft_selectdata(data, 'rpt', cfg.trials);
  refdata = ft_selectdata(refdata, 'rpt', cfg.trials);
end

refchan = ft_channelselection(cfg.refchannel, refdata.label);
refindx = match_str(refdata.label, refchan);
megchan = ft_channelselection(cfg.channel, data.label);
megindx = match_str(data.label, megchan);
if ~computeweights
  %check whether the weight table contains the specified references
  %ensure the ordering of the meg-data to be consistent with the weights
  %ensure the ordering of the ref-data to be consistent with the weights
  
  [i1,i2] = match_str(refchan, cfg.pca.reflabel);
  [i3,i4] = match_str(megchan, cfg.pca.label);
  if length(i2)~=length(cfg.pca.reflabel),
    error('you specified fewer references to use as there are in the precomputed weight table');
  end
  
  refindx = refindx(i1);
  megindx = megindx(i3);
  cfg.pca.w = cfg.pca.w(i4,i2);
  cfg.pca.label   = cfg.pca.label(i4);
  cfg.pca.reflabel= cfg.pca.reflabel(i2);
  if isfield(cfg.pca, 'rotmat'),
    cfg.pca = rmfield(cfg.pca, 'rotmat'); %dont know
  end
else
  %do nothing
end
nmeg    = length(megindx);
nref    = length(refindx);

if ischar(cfg.truncate) && strcmp(cfg.truncate, 'no'),
  cfg.truncate = length(refindx);
elseif ischar(cfg.truncate) || (cfg.truncate>1 && cfg.truncate/round(cfg.truncate)~=1) || ...
       cfg.truncate>length(refindx),
  error('cfg.truncate should be either ''no'', an integer number <= the number of references, or a number between 0 and 1');
  %FIXME the default truncation applied by 4D is 1x10^-8
end

ntrl   = length(data.trial);
nchan  = length(data.label);
nsmp   = cellfun('size', data.trial, 2);

%compute and remove mean from data
if strcmp(cfg.blc, 'yes'),
  m             = cellmean(data.trial,    2);
  data.trial    = cellvecadd(data.trial, -m);
  m             = cellmean(refdata.trial,    2);
  refdata.trial = cellvecadd(refdata.trial, -m);
end

%compute std of data before
stdpre = cellstd(data.trial, 2);

if computeweights,
  %zscore
  if strcmp(cfg.zscore, 'yes'),
    [refdata.trial, sdref] = cellzscore(refdata.trial, 2, 0); %forced demeaned already
  else
    sdref = ones(nref,1);
  end
  
  %compute covariance of refchannels and do svd
  crefdat = cellcov(refdata.trial, [], 2, 0);
  [u,s,v] = svd(crefdat);
  
  %determine the truncation and rotation
  if cfg.truncate<1
    %keep all singular vectors with singular values >= cfg.truncate*s(1,1)
    s1   = s./max(s(:));
    keep = find(diag(s1)>cfg.truncate);
  else
    keep = 1:cfg.truncate;
  end
  rotmat = u(:, keep)';
  
  %rotate the refdata
  refdata.trial = cellfun(@mtimes, repmat({rotmat}, 1, ntrl), refdata.trial, 'UniformOutput', 0);
  
  %project megdata onto the orthogonal basis
  nom   = cellcov(data.trial, refdata.trial, 2, 0);
  denom = cellcov(refdata.trial, [],         2, 0);
  rw    = (pinv(denom)*nom')';
  
  %subtract projected data
  for k = 1:ntrl
    data.trial{k} = data.trial{k} - rw*refdata.trial{k};
  end
  
  %rotate back and 'unscale'
  pca.w  = rw*rotmat*diag(1./sdref);
  pca.label    = data.label;
  pca.reflabel = refdata.label;
  pca.rotmat   = rotmat;
  cfg.pca      = pca;
else
  fprintf('applying precomputed weights to the data\n');
  pca = cfg.pca;
  for k = 1:ntrl
    data.trial{k} = data.trial{k} - pca.w*refdata.trial{k};
  end
end

%compute std of data after
stdpst = cellstd(data.trial, 2);

%blc
if strcmp(cfg.blc, 'yes'),
  m          = cellmean(data.trial, 2);
  data.trial = cellvecadd(data.trial, -m);
end

%apply weights to the gradiometer-array
if isfield(data, 'grad')
  montage     = [];
  labelnew    = pca.label;
  nlabelnew   = length(labelnew);
  
  %add columns of refchannels not yet present in labelnew
  %[id, i1]  = setdiff(pca.reflabel, labelnew);
  %labelorg  = [labelnew; pca.reflabel(sort(i1))];
  labelorg  = data.grad.label;
  nlabelorg = length(labelorg);

  %start with identity
  montage.tra = eye(nlabelorg);
  
  %subtract weights
  [i1, i2]  = match_str(labelorg, pca.reflabel);
  [i3, i4]  = match_str(labelorg, pca.label);
  montage.tra(i3,i1) = montage.tra(i3,i1) - pca.w(i4,i2);
  montage.labelorg  = labelorg;
  montage.labelnew  = labelorg;
  data.grad         = ft_apply_montage(data.grad, montage, 'keepunused', 'yes');
  
  if isfield(data.grad, 'balance'),
    mnt = fieldnames(data.grad.balance);
    sel = strmatch('pca', mnt);
    if isempty(sel),
      sel = zeros(0,1);
    end  
    name    = ['pca',num2str(length(sel)+1)];
    if ~strcmp(data.grad.balance.current, 'none'),
      current = getsubfield(data.grad.balance, data.grad.balance.current);
      data.grad.balance.current = [name,'_',data.grad.balance.current];
    else
      current.tra      = eye(size(data.grad.tra,1));
      current.labelorg = data.grad.label;
      current.labelnew = data.grad.label;
      data.grad.balance.current = name;
    end
    newcurrent = ft_apply_montage(current, montage, 'keepunused', 'yes');
    data.grad.balance = setsubfield(data.grad.balance, data.grad.balance.current, newcurrent);
  else
    balance = struct('pca1', montage);
    data.grad.balance = balance;
    data.grad.balance.current = 'pca1';
  end
else
  warning('weights have been applied to the data only, not to the sensors');
end

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% add the version details of this function call to the configuration
cfg.version.name = mfilename('fullpath');
cfg.version.id   = '$Id: ft_denoise_pca.m,v 1.9 2009/10/12 14:50:59 jansch Exp $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output
data.cfg = cfg;

%-----cellcov
function [c] = cellcov(x, y, dim, flag)

% [C] = CELLCOV(X, DIM) computes the covariance, across all cells in x along 
% the dimension dim. When there are three inputs, covariance is computed between
% all cells in x and y
% 
% X (and Y) should be linear cell-array(s) of matrices for which the size in at 
% least one of the dimensions should be the same for all cells 

if nargin==2,
  flag = 1;
  dim  = y;
  y    = [];
elseif nargin==3,
  flag = 1;
end

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellmean');
end

if nargin==1,
  scx1 = cellfun('size', x, 1);
  scx2 = cellfun('size', x, 2);
  if     all(scx2==scx2(1)), dim = 2; %let second dimension prevail
  elseif all(scx1==scx1(1)), dim = 1;
  else   error('no dimension to compute covariance for');
  end
end

if flag,
  mx   = cellmean(x, 2);
  x    = cellvecadd(x, -mx);
  if ~isempty(y),
    my = cellmean(y, 2);
    y  = cellvecadd(y, -my);
  end
end

nx   = max(nx);
nsmp = cellfun('size', x, dim);
if isempty(y), 
  csmp = cellfun(@covc, x, repmat({dim},1,nx), 'UniformOutput', 0);
else
  csmp = cellfun(@covc, x, y, repmat({dim},1,nx), 'UniformOutput', 0);
end
nc   = size(csmp{1});
c    = sum(reshape(cell2mat(csmp), [nc(1) nc(2) nx]), 3)./sum(nsmp); 

function [c] = covc(x, y, dim)

if nargin==2,
  dim = y;
  y   = x;
end

if dim==1,
  c = x'*y;
elseif dim==2,
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
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellmean');
end

if nargin==1,
  scx1 = cellfun('size', x, 1);
  scx2 = cellfun('size', x, 2);
  if     all(scx2==scx2(1)), dim = 2; %let second dimension prevail
  elseif all(scx1==scx1(1)), dim = 1;
  else   error('no dimension to compute mean for');
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
% be subtracted first (default behaviour, but to save time on already demeaned data, it
% can be set to 0).

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellstd');
end

if nargin<2,
  scx1 = cellfun('size', x, 1);
  scx2 = cellfun('size', x, 2);
  if     all(scx2==scx2(1)), dim = 2; %let second dimension prevail
  elseif all(scx1==scx1(1)), dim = 1;
  else   error('no dimension to compute mean for');
  end
elseif nargin==2,
  flag = 1;
end

if flag,
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
if isempty(bsxfun_exists); 
    bsxfun_exists=(exist('bsxfun')==5); 
    if ~bsxfun_exists; 
        error('bsxfun not found.');
    end
end

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellmean');
end

if ~iscell(v),
  v = repmat({v}, nx);
end

sx1 = cellfun('size', x, 1);
sx2 = cellfun('size', x, 2);
sv1 = cellfun('size', v, 1);
sv2 = cellfun('size', v, 2);
if all(sx1==sv1) && all(sv2==1),    
  dim = mat2cell([ones(length(sx2),1) sx2(:)]', repmat(2,nx(1),1), repmat(1,nx(2),1)); 
elseif all(sx2==sv2) && all(sv1==1),
  dim = mat2cell([sx1(:) ones(length(sx1),1)]', repmat(2,nx(1),1), repmat(1,nx(2),1));
elseif all(sv1==1) && all(sv2==1),
  dim = mat2cell([sx1(:) sx2(:)]'', nx(1), nx(2));
else   error('inconsistent input');
end  

y  = cellfun(@bsxfun, repmat({@plus}, nx), x, v, 'UniformOutput', 0);
%y = cellfun(@vplus, x, v, dim, 'UniformOutput', 0);

function y = vplus(x, v, dim)

y = x + repmat(v, dim);

%-----cellvecmult
function [y] = cellvecmult(x, v)

% [Y]= CELLVECMULT(X, V) - multiply vectors in cell-array V
% to all rows or columns of each matrix in cell-array X
% V can be a vector or a cell-array of vectors

% check once and for all to save time
persistent bsxfun_exists;
if isempty(bsxfun_exists); 
    bsxfun_exists=(exist('bsxfun')==5); 
    if ~bsxfun_exists; 
        error('bsxfun not found.');
    end
end

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellmean');
end

if ~iscell(v),
  v = repmat({v}, nx);
end

sx1 = cellfun('size', x, 1);
sx2 = cellfun('size', x, 2);
sv1 = cellfun('size', v, 1);
sv2 = cellfun('size', v, 2);
if all(sx1==sv1) && all(sv2==1),    
elseif all(sx2==sv2) && all(sv1==1),
elseif all(sv1==1) && all(sv2==1),
else   error('inconsistent input');
end  

y  = cellfun(@bsxfun, repmat({@times}, nx), x, v, 'UniformOutput', 0);

%-----cellzscore
function [z, sd, m] = cellzscore(x, dim, flag)

% [Z, SD] = CELLZSCORE(X, DIM, FLAG) computes the zscore, across all cells in x along 
% the dimension dim, normalising by the total number of samples 
% 
% X should be an linear cell-array of matrices for which the size in at 
% least one of the dimensions should be the same for all cells. If flag==1, the mean will
% be subtracted first (default behaviour, but to save time on already demeaned data, it
% can be set to 0). SD is a vector containing the standard deviations, used for the normalisation.

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellstd');
end

if nargin<2,
  scx1 = cellfun('size', x, 1);
  scx2 = cellfun('size', x, 2);
  if     all(scx2==scx2(1)), dim = 2; %let second dimension prevail
  elseif all(scx1==scx1(1)), dim = 1;
  else   error('no dimension to compute mean for');
  end
elseif nargin==2,
  flag = 1;
end

if flag,
  m    = cellmean(x, dim);
  x    = cellvecadd(x, -m);
end

sd   = cellstd(x, dim, 0);
z    = cellvecmult(x, 1./sd);

