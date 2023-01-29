function [parcel] = ft_sourceparcellate(cfg, source, parcellation)

% FT_SOURCEPARCELLATE combines the source-reconstruction parameters over the parcels, for
% example by averaging all the values in the anatomically or functionally labeled parcel.
%
% Use as
%    output = ft_sourceparcellate(cfg, source, parcellation)
% where the input source is a 2D surface-based or 3-D voxel-based source grid that was for
% example obtained from FT_SOURCEANALYSIS or FT_COMPUTE_LEADFIELD. The input parcellation is
% described in detail in FT_DATATYPE_PARCELLATION (2-D) or FT_DATATYPE_SEGMENTATION (3-D) and
% can be obtained from FT_READ_ATLAS or from a custom parcellation/segmentation for your
% individual subject. The output is a channel-based representation with the combined (e.g.
% averaged) representation of the source parameters per parcel.
%
% The configuration "cfg" is a structure that can contain the following fields
%   cfg.method       = string, method to combine the values, see below (default = 'mean')
%   cfg.parcellation = string, fieldname that contains the desired parcellation
%   cfg.parameter    = cell-array with strings, fields that should be parcellated (default = 'all')
%
% The values within a parcel or parcel-combination can be combined with different methods:
%   'mean'      compute the mean
%   'median'    compute the median (unsupported for fields that are represented in a cell-array)
%   'eig'       compute the largest eigenvector
%   'min'       take the minimal value
%   'max'       take the maximal value
%   'maxabs'    take the signed maxabs value
%   'std'       take the standard deviation
%
% See also FT_SOURCEANALYSIS, FT_DATATYPE_PARCELLATION, FT_DATATYPE_SEGMENTATION

% Copyright (C) 2012-2021, Robert Oostenveld
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
ft_preamble loadvar source parcellation
ft_preamble provenance source parcellation

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% get the defaults
cfg.parcellation = ft_getopt(cfg, 'parcellation');
cfg.parameter    = ft_getopt(cfg, 'parameter', 'all');
cfg.method       = ft_getopt(cfg, 'method', 'mean'); % can be mean, min, max, svd
cfg.feedback     = ft_getopt(cfg, 'feedback', 'text');

% the data can be passed as input argument or can be read from disk
hasparcellation = exist('parcellation', 'var');

if ischar(cfg.parameter)
  cfg.parameter = {cfg.parameter};
end

if hasparcellation
  % the parcellation is specified as separate structure
else
  % the parcellation is represented in the source structure itself
  parcellation = source;
end

% keep the transformation matrix
if isfield(parcellation, 'transform')
  transform = parcellation.transform;
else
  transform = [];
end

% ensure it is a parcellation, not a segmentation
parcellation = ft_checkdata(parcellation, 'datatype', 'parcellation', 'parcellationstyle', 'indexed', 'hasunit', 'yes');

% keep the transformation matrix
if ~isempty(transform)
  parcellation.transform = transform;
end

% ensure it is a source, not a volume
source = ft_checkdata(source, 'datatype', 'source', 'insidestyle', 'logical', 'hasunit', 'yes');

% ensure that the source and the parcellation are anatomically consistent
if ~strcmp(source.unit, parcellation.unit)
  ft_error('the units of the source and parcellation structure are not consistent, please use FT_SOURCEINTERPOLATE');
end

% ensure that the source and the parcellation are anatomically consistent
tolerance = 0.1 * ft_scalingfactor('mm', source.unit);
if ~isalmostequal(source.pos, parcellation.pos, 'abstol', tolerance)
  ft_error('the positions of the source and parcellation structure are not consistent, please use FT_SOURCEINTERPOLATE');
end

if isempty(cfg.parcellation)
  % determine the first field that can be used for the parcellation
  fn = fieldnames(parcellation);
  for i=1:numel(fn)
    if isfield(parcellation, [fn{i} 'label'])
      ft_warning('using "%s" for the parcellation', fn{i});
      cfg.parcellation = fn{i};
      break
    end
  end
end

if isempty(cfg.parcellation)
  ft_error('you should specify the field containing the parcellation');
end

% determine the fields and corresponding dimords to work on
fn = fieldnames(source);
fn = setdiff(fn, {'pos', 'tri', 'inside', 'outside', 'time', 'freq', 'dim', 'transform', 'unit', 'coordsys', 'cfg', 'hdr'}); % remove fields that do not represent the data
fn = fn(cellfun(@isempty, regexp(fn, 'dimord'))); % remove dimord fields
fn = fn(cellfun(@isempty, regexp(fn, 'label'))); % remove label fields
dimord = cell(size(fn));
for i=1:numel(fn)
  dimord{i} = getdimord(source, fn{i});
end

if any(strcmp(cfg.parameter, 'all'))
  cfg.parameter = fn;
else
  [inside, i1, i2] = intersect(cfg.parameter, fn);
  [outside       ] = setdiff(cfg.parameter, fn);
  if ~isempty(outside)
    ft_warning('\nparameter "%s" cannot be parcellated', outside{:});
  end
  cfg.parameter = fn(i2);
  fn     = fn(i2);
  dimord = dimord(i2);
end

% although it is technically feasible, don't parcellate the parcellation itself
sel    = ~strcmp(cfg.parcellation, fn);
fn     = fn(sel);
dimord = dimord(sel);

if numel(fn)==0
  ft_error('there are no source parameters that can be parcellated');
end

% get the parcellation and the labels that go with it
tissue      = parcellation.(cfg.parcellation);
tissuelabel = parcellation.([cfg.parcellation 'label']);
ntissue     = length(tissuelabel);

if isfield(source, 'inside')
  % determine the conjunction of the parcellation and the inside source points
  n0 = numel(source.inside);
  n1 = sum(source.inside(:));
  n2 = sum(tissue(:)~=0);
  fprintf('there are in total %d positions, %d positions are inside the brain, %d positions have a label\n', n0, n1, n2);
  fprintf('%d of the positions inside the brain have a label\n',        sum(tissue(source.inside)~=0));
  fprintf('%d of the labeled positions are inside the brain\n',         sum(source.inside(tissue(:)~=0)));
  fprintf('%d of the positions inside the brain do not have a label\n', sum(tissue(source.inside)==0));
  % discard the positions outside the brain and the positions in the brain that do not have a label
  tissue(~source.inside) = 0;
end

% start preparing the output data structure
parcel       = keepfields(source, {'freq','time','cumtapcnt'});
parcel.label = tissuelabel;

for i=1:numel(fn)
  % parcellate each of the desired parameters
  dat = source.(fn{i});
  siz = getdimsiz(source, fn{i});
  siz(contains(tokenize(dimord{i},'_'),'pos')) = ntissue;
  
  if startsWith(dimord{i}, '{pos_pos}')
    fprintf('creating %d*%d parcel combinations for parameter %s by taking the %s\n', numel(tissuelabel), numel(tissuelabel), fn{i}, cfg.method);
    tmp = zeros(siz);
    ft_progress('init', cfg.feedback, 'computing parcellation');
    k = 0;
    K = numel(tissuelabel)^2;
    for j1=1:numel(tissuelabel)
      for j2=1:numel(tissuelabel)
        k = k + 1;
        ft_progress(k/K, 'computing parcellation for %s combined with %s', tissuelabel{j1}, tissuelabel{j2});
        switch cfg.method
          case 'mean'
            tmp(j1,j2,:,:) = cellmean(dat(tissue==j1,tissue==j2));
          case 'median'
            tmp(j1,j2,:,:) = cellmedian(dat(tissue==j1,tissue==j2));
          case 'min'
            tmp(j1,j2,:,:) = cellmin(dat(tissue==j1,tissue==j2));
          case 'max'
            tmp(j1,j2,:,:) = cellmax(dat(tissue==j1,tissue==j2));
          case 'eig'
            tmp(j1,j2,:,:) = celleig(dat(tissue==j1,tissue==j2));
          case 'std'
            tmp(j1,j2,:,:) = cellstd(dat(tissue==j1,tissue==j2));
          otherwise
            ft_error('method %s not implemented for %s', cfg.method, dimord{i});
        end % switch
      end % for j2
    end % for j1
    ft_progress('close');
    
  elseif startsWith(dimord{i}, '{pos}')
    fprintf('creating %d parcels for parameter %s by taking the %s\n', numel(tissuelabel), fn{i}, cfg.method);
    tmp = zeros(siz);
    ft_progress('init', cfg.feedback, 'computing parcellation');
    for j=1:numel(tissuelabel)
      ft_progress(j/numel(tissuelabel), 'computing parcellation for %s', tissuelabel{j});
      switch cfg.method
        case 'mean'
          tmp(j,:,:) = cellmean(dat(tissue==j));
        case 'median'
          tmp(j,:,:) = cellmedian(dat(tissue==j));
        case 'min'
          tmp(j,:,:) = cellmin(dat(tissue==j));
        case 'max'
          tmp(j,:,:) = cellmax(dat(tissue==j));
        case 'eig'
          tmp(j,:,:) = celleig(dat(tissue==j));
        case 'std'
          tmp(j,:,:) = cellstd(dat(tissue==j));
        otherwise
          ft_error('method %s not implemented for %s', cfg.method, dimord{i});
      end % switch
    end % for
    ft_progress('close');
    
  elseif startsWith(dimord{i}, 'pos_pos')
    fprintf('creating %d*%d parcel combinations for parameter %s by taking the %s\n', numel(tissuelabel), numel(tissuelabel), fn{i}, cfg.method);
    siz     = size(dat);
    siz(1)  = ntissue;
    siz(2)  = ntissue;
    tmp     = nan(siz);
    ft_progress('init', cfg.feedback, 'computing parcellation');
    k = 0;
    K = numel(tissuelabel)^2;
    for j1=1:numel(tissuelabel)
      for j2=1:numel(tissuelabel)
        k = k + 1;
        ft_progress(k/K, 'computing parcellation for %s combined with %s', tissuelabel{j1}, tissuelabel{j2});
        switch cfg.method
          case 'mean'
            tmp(j1,j2,:) = arraymean2(dat(tissue==j1,tissue==j2,:));
          case 'median'
            tmp(j1,j2,:) = arraymedian2(dat(tissue==j1,tissue==j2,:));
          case 'min'
            tmp(j1,j2,:) = arraymin2(dat(tissue==j1,tissue==j2,:));
          case 'max'
            tmp(j1,j2,:) = arraymax2(dat(tissue==j1,tissue==j2,:));
          case 'eig'
            tmp(j1,j2,:) = arrayeig2(dat(tissue==j1,tissue==j2,:));
          case 'maxabs'
            tmp(j1,j2,:) = arraymaxabs2(dat(tissue==j1,tissue==j2,:));
          case 'std'
            tmp(j1,j2,:) = arraystd2(dat(tissue==j1,tissue==j2,:));
          otherwise
            ft_error('method %s not implemented for %s', cfg.method, dimord{i});
        end % switch
      end % for j2
    end % for j1
    ft_progress('close');
    
  elseif startsWith(dimord{i}, 'pos')
    fprintf('creating %d parcels for %s by taking the %s\n', numel(tissuelabel), fn{i}, cfg.method);
    siz     = size(dat);
    siz(1)  = ntissue;
    tmp     = nan(siz);
    ft_progress('init', cfg.feedback, 'computing parcellation');
    for j=1:numel(tissuelabel)
      ft_progress(j/numel(tissuelabel), 'computing parcellation for %s', tissuelabel{j});
      switch cfg.method
        case 'mean'
          tmp(j,:) = arraymean1(dat(tissue==j,:));
        case 'mean_thresholded'
          cfg.mean = ft_getopt(cfg, 'mean', struct('threshold', []));
          if isempty(cfg.mean.threshold)
            ft_error('when cfg.method = ''mean_thresholded'', you should specify a cfg.mean.threshold');
          end
          if numel(cfg.mean.threshold)==size(dat,1)
            % assume one threshold per vertex
            threshold = cfg.mean.threshold(tissue==j,:);
          else
            threshold = cfg.mean.threshold;
          end
          tmp(j,:) = arraymean1(dat(tissue==j,:), threshold);
        case 'median'
          tmp(j,:) = arraymedian1(dat(tissue==j,:));
        case 'min'
          tmp(j,:) = arraymin1(dat(tissue==j,:));
        case 'max'
          tmp(j,:) = arraymax1(dat(tissue==j,:));
        case 'maxabs'
          tmp(j,:) = arraymaxabs1(dat(tissue==j,:));
        case 'eig'
          tmp(j,:) = arrayeig1(dat(tissue==j,:));
        case 'std'
          tmp(j,:)  = arraystd1(dat(tissue==j,:));
        otherwise
          ft_error('method %s not implemented for %s', cfg.method, dimord{i});
      end % switch
    end % for
    ft_progress('close');
    
  else
    ft_error('unsupported dimord %s', dimord{i})
    
  end % if pos, pos_pos, {pos}, etc.
  
  % update the dimord, use chan rather than pos
  % this makes it look just like timelock or freq data
  tok = tokenize(dimord{i}, '_');
  tok(strcmp(tok,  'pos' )) = {'chan'}; % replace pos by chan
  tok(strcmp(tok, '{pos}')) = {'chan'}; % replace pos by chan
  tok(strcmp(tok, '{pos'))  = {'chan'}; % replace pos by chan
  tok(strcmp(tok, 'pos}'))  = {'chan'}; % replace pos by chan
  
  % squeeze out any singleton oris
  siz  = [size(tmp) 1]; % add trailing singleton to be sure
  oris = contains(tok, 'ori') & siz(1:numel(tok))==1;
  siz(oris) = [];
  tmp = reshape(tmp, siz);
  tok(oris) = [];
  
  tmpdimord = sprintf('%s_', tok{:});
  tmpdimord = tmpdimord(1:end-1);         % exclude the last _
  
  % store the results in the output structure
  parcel.(fn{i})            = tmp;
  parcel.([fn{i} 'dimord']) = tmpdimord;
  
  % to avoid confusion
  clear dat tmp tmpdimord j j1 j2
end % for each of the fields that should be parcellated

% a brainordinate is a brain location that is specified by either a surface vertex (node) or a volume voxel
parcel.brainordinate = keepfields(parcellation, {'pos', 'tri', 'dim', 'transform'}); % keep the information about the geometry
fn = fieldnames(parcellation);
for i=1:numel(fn)
  if isfield(parcellation, [fn{i} 'label'])
    % keep each of the labeled fields from the parcellation
    parcel.brainordinate.( fn{i}         ) = parcellation.( fn{i}         );
    parcel.brainordinate.([fn{i} 'label']) = parcellation.([fn{i} 'label']);
  end
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous   source parcellation
ft_postamble provenance parcel
ft_postamble history    parcel
ft_postamble savevar    parcel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS to compute something over the first dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = arraymean1(x, threshold)

if nargin==1
  y = mean(x,1);
else
  if numel(threshold)==1
    % scalar comparison is possible
  elseif size(threshold,1) == size(x,1)
    % assume threshold to be column vector
    threshold = repmat(threshold, [1, size(x,2)]);
  end
  sel = sum(x>threshold,2);
  if ~isempty(sel)
    y = mean(x(sel>0,:),1);
  else
    y = nan+zeros(1,size(x,2));
  end
end

function y = arraymedian1(x)
if ~isempty(x)
  y = median(x,1);
else
  y = nan(1,size(x,2));
end

function y = arraymin1(x)
if ~isempty(x)
  y = min(x,[], 1);
else
  y = nan(1,size(x,2));
end

function y = arraymax1(x)
if ~isempty(x)
  y = max(x,[], 1);
else
  y = nan(1,size(x,2));
end

function y = arrayeig1(x)
if ~isempty(x)
  siz = size(x);
  x = reshape(x, siz(1), prod(siz(2:end)));
  [u, s, v] = svds(x, 1);         % x = u * s * v'
  y = s(1,1) * v(:,1);            % retain the largest eigenvector with appropriate scaling
  y = reshape(y, [siz(2:end) 1]); % size should have at least two elements
else
  siz = size(x);
  y   = nan([siz(2:end) 1]);
end

function y = arraymaxabs1(x)
if ~isempty(x)
  % take the value that is at max(abs(x))
  [dum,ix] = max(abs(x), [], 1);
  y        = x(ix);
else
  y = nan(1,size(x,2));
end

function y = arraystd1(x)
if ~isempty(x)
  y = std(x,0, 1);
else
  y = nan(1,size(x, 2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS to compute something over the first two dimensions
% all of these functions should be implemented the same
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = arraymean2(x)
siz = size(x);
x = reshape(x, [siz(1)*siz(2) siz(3:end) 1]); % simplify it into a single dimension
y = arraymean1(x);

function y = arraymedian2(x)
siz = size(x);
x = reshape(x, [siz(1)*siz(2) siz(3:end) 1]); % simplify it into a single dimension
y = arraymedian1(x);

function y = arraymin2(x)
siz = size(x);
x = reshape(x, [siz(1)*siz(2) siz(3:end) 1]); % simplify it into a single dimension
y = arraymin1(x);

function y = arraymax2(x)
siz = size(x);
x = reshape(x, [siz(1)*siz(2) siz(3:end) 1]); % simplify it into a single dimension
y = arraymax1(x);

function y = arrayeig2(x)
siz = size(x);
x = reshape(x, [siz(1)*siz(2) siz(3:end) 1]); % simplify it into a single dimension
y = arrayeig1(x);

function y = arraymaxabs2(x)
siz = size(x);
x = reshape(x, [siz(1)*siz(2) siz(3:end) 1]); % simplify it into a single dimension
y = arraymaxabs1(x);

function y = arraystd2(x)
siz = size(x);
x = reshape(x, [siz(1)*siz(2) siz(3:end) 1]); % simplify it into a single dimension
y = arraystd1(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS for doing something over all elements of a cell-array
% add a singleton dimension, concatenate into an array, and do the computatioon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = cellmean(x)
siz = size(x{1});
for i=1:numel(x)
  x{i} = reshape(x{i}, [1 siz]);
end
x = cat(1, x{:});
y = arraymean1(x);

function y = cellmedian(x)
siz = size(x{1});
for i=1:numel(x)
  x{i} = reshape(x{i}, [1 siz]);
end
x = cat(1, x{:});
y = arraymedian1(x);

function y = cellmin(x)
siz = size(x{1});
for i=1:numel(x)
  x{i} = reshape(x{i}, [1 siz]);
end
x = cat(1, x{:});
y = arraymin1(x);

function y = cellmax(x)
siz = size(x{1});
for i=1:numel(x)
  x{i} = reshape(x{i}, [1 siz]);
end
x = cat(1, x{:});
y = arraymax1(x);

function y = celleig(x)
siz = size(x{1});
for i=1:numel(x)
  x{i} = reshape(x{i}, [1 siz]);
end
x = cat(1, x{:});
y = arrayeig1(x);

function y = cellstd(x)
siz = size(x{1});
for i=1:numel(x)
  x{i} = reshape(x{i}, [1 siz]);
end
x = cat(1, x{:});
y = arraystd1(x);
