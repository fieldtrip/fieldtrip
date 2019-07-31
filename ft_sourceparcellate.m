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

% Copyright (C) 2012-2013, Robert Oostenveld
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
ft_preamble trackconfig

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
parcellation = ft_checkdata(parcellation, 'datatype', 'parcellation', 'parcellationstyle', 'indexed');

% keep the transformation matrix
if ~isempty(transform)
  parcellation.transform = transform;
end

% ensure it is a source, not a volume
source = ft_checkdata(source, 'datatype', 'source', 'inside', 'logical');

% ensure that the source and the parcellation are anatomically consistent
if ~isalmostequal(source.pos, parcellation.pos, 'abstol', 1000000*eps)
  ft_error('the source positions are not consistent with the parcellation, please use FT_SOURCEINTERPOLATE');
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
seg      = parcellation.(cfg.parcellation);
seglabel = parcellation.([cfg.parcellation 'label']);
nseg     = length(seglabel);

if isfield(source, 'inside')
  % determine the conjunction of the parcellation and the inside source points
  n0 = numel(source.inside);
  n1 = sum(source.inside(:));
  n2 = sum(seg(:)~=0);
  fprintf('there are in total %d positions, %d positions are inside the brain, %d positions have a label\n', n0, n1, n2);
  fprintf('%d of the positions inside the brain have a label\n',        sum(seg(source.inside)~=0));
  fprintf('%d of the labeled positions are inside the brain\n',         sum(source.inside(seg(:)~=0)));
  fprintf('%d of the positions inside the brain do not have a label\n', sum(seg(source.inside)==0));
  % discard the positions outside the brain and the positions in the brain that do not have a label
  seg(~source.inside) = 0;
end

% start preparing the output data structure
parcel       = keepfields(source, {'freq','time','cumtapcnt'});
parcel.label = seglabel;

for i=1:numel(fn)
  % parcellate each of the desired parameters
  dat = source.(fn{i});
  
  if strncmp('{pos_pos}', dimord{i}, 9)
    fprintf('creating %d*%d parcel combinations for parameter %s by taking the %s\n', numel(seglabel), numel(seglabel), fn{i}, cfg.method);
    tmp = cell(nseg, nseg);
    ft_progress('init', cfg.feedback, 'computing parcellation');
    k = 0;
    K = numel(seglabel)^2;
    for j1=1:numel(seglabel)
      for j2=1:numel(seglabel)
        k = k + 1;
        ft_progress(k/K, 'computing parcellation for %s combined with %s', seglabel{j1}, seglabel{j2});
        switch cfg.method
          case 'mean'
            tmp{j1,j2} = cellmean2(dat(seg==j1,seg==j2,:));
          case 'median'
            tmp{j1,j2} = cellmedian2(dat(seg==j1,seg==j2,:));
          case 'min'
            tmp{j1,j2} = cellmin2(dat(seg==j1,seg==j2,:));
          case 'max'
            tmp{j1,j2} = cellmax2(dat(seg==j1,seg==j2,:));
          case 'eig'
            tmp{j1,j2} = celleig2(dat(seg==j1,seg==j2,:));
          case 'std'
            tmp{j1,j2} = cellstd2(dat(seg==j1,seg==j2,:));
          otherwise
            ft_error('method %s not implemented for %s', cfg.method, dimord{i});
        end % switch
      end % for j2
    end % for j1
    ft_progress('close');
    
  elseif strncmp('{pos}', dimord{i}, 5)
    fprintf('creating %d parcels for parameter %s by taking the %s\n', numel(seglabel), fn{i}, cfg.method);
    tmp = cell(nseg, 1);
    ft_progress('init', cfg.feedback, 'computing parcellation');
    for j=1:numel(seglabel)
      ft_progress(j/numel(seglabel), 'computing parcellation for %s', seglabel{j});
      switch cfg.method
        case 'mean'
          tmp{j} = cellmean1(dat(seg==j));
        case 'median'
          tmp{j} = cellmedian1(dat(seg==j));
        case 'min'
          tmp{j} = cellmin1(dat(seg==j));
        case 'max'
          tmp{j} = cellmax1(dat(seg==j));
        case 'eig'
          tmp{j} = celleig1(dat(seg==j));
        case 'std'
          tmp{j} = cellstd1(dat(seg==j));
        otherwise
          ft_error('method %s not implemented for %s', cfg.method, dimord{i});
      end % switch
    end % for
    ft_progress('close');
    
  elseif strncmp('pos_pos', dimord{i}, 7)
    fprintf('creating %d*%d parcel combinations for parameter %s by taking the %s\n', numel(seglabel), numel(seglabel), fn{i}, cfg.method);
    siz     = size(dat);
    siz(1)  = nseg;
    siz(2)  = nseg;
    tmp     = nan(siz);
    ft_progress('init', cfg.feedback, 'computing parcellation');
    k = 0;
    K = numel(seglabel)^2;
    for j1=1:numel(seglabel)
      for j2=1:numel(seglabel)
        k = k + 1;
        ft_progress(k/K, 'computing parcellation for %s combined with %s', seglabel{j1}, seglabel{j2});
        switch cfg.method
          case 'mean'
            tmp(j1,j2,:) = arraymean2(dat(seg==j1,seg==j2,:));
          case 'median'
            tmp(j1,j2,:) = arraymedian2(dat(seg==j1,seg==j2,:));
          case 'min'
            tmp(j1,j2,:) = arraymin2(dat(seg==j1,seg==j2,:));
          case 'max'
            tmp(j1,j2,:) = arraymax2(dat(seg==j1,seg==j2,:));
          case 'eig'
            tmp(j1,j2,:) = arrayeig2(dat(seg==j1,seg==j2,:));
          case 'maxabs'
            tmp(j1,j2,:) = arraymaxabs2(dat(seg==j1,seg==j2,:));
          case 'std'
            tmp(j1,j2,:) = arraystd2(dat(seg==j1,seg==j2,:));
          otherwise
            ft_error('method %s not implemented for %s', cfg.method, dimord{i});
        end % switch
      end % for j2
    end % for j1
    ft_progress('close');
    
  elseif strncmp('pos', dimord{i}, 3)
    fprintf('creating %d parcels for %s by taking the %s\n', numel(seglabel), fn{i}, cfg.method);
    siz     = size(dat);
    siz(1)  = nseg;
    tmp     = nan(siz);
    ft_progress('init', cfg.feedback, 'computing parcellation');
    for j=1:numel(seglabel)
      ft_progress(j/numel(seglabel), 'computing parcellation for %s', seglabel{j});
      switch cfg.method
        case 'mean'
          tmp(j,:) = arraymean1(dat(seg==j,:));
        case 'mean_thresholded'
          cfg.mean = ft_getopt(cfg, 'mean', struct('threshold', []));
          if isempty(cfg.mean.threshold)
            ft_error('when cfg.method = ''mean_thresholded'', you should specify a cfg.mean.threshold');
          end
          if numel(cfg.mean.threshold)==size(dat,1)
            % assume one threshold per vertex
            threshold = cfg.mean.threshold(seg==j,:);
          else
            threshold = cfg.mean.threshold;
          end
          tmp(j,:) = arraymean1(dat(seg==j,:), threshold);
        case 'median'
          tmp(j,:) = arraymedian1(dat(seg==j,:));
        case 'min'
          tmp(j,:) = arraymin1(dat(seg==j,:));
        case 'max'
          tmp(j,:) = arraymax1(dat(seg==j,:));
        case 'maxabs'
          tmp(j,:) = arraymaxabs1(dat(seg==j,:));
        case 'eig'
          tmp(j,:) = arrayeig1(dat(seg==j,:));
        case 'std'
          tmp(j,:)  = arraystd1(dat(seg==j,:));
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
  tok(strcmp(tok,  'pos' )) = { 'chan' }; % replace pos by chan
  tok(strcmp(tok, '{pos}')) = {'{chan}'}; % replace pos by chan
  tok(strcmp(tok, '{pos'))  = {'{chan' }; % replace pos by chan
  tok(strcmp(tok, 'pos}'))  = { 'chan}'}; % replace pos by chan
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

ft_postamble debug
ft_postamble trackconfig
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
    y   = mean(x(sel>0,:),1);
  else
    y   = nan+zeros(1,size(x,2));
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
% take the value that is at max(abs(x))
siz = size(x);
x = reshape(x, [siz(1)*siz(2) siz(3:end) 1]); % simplify it into a single dimension
[dum,ix] = max(abs(x), [], 1);
y        = x(ix);

function y = arraystd2(x)
siz = size(x);
x = reshape(x, [siz(1)*siz(2) siz(3:end) 1]); % simplify it into a single dimension
y = arraystd1(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS for doing something over the first dimension of a cell-array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = cellmean1(x)
siz = size(x);
if siz(1)==1 && siz(2)>1
  siz([2 1]) = siz([1 2]);
  x = reshape(x, siz);
end
y = x{1};
n = 1;
for i=2:siz(1)
  y = y + x{i};
  n = n + 1;
end
y = y/n;

function y = cellmedian1(x)
siz = size(x);
if siz(1)==1 && siz(2)>1
  siz([2 1]) = siz([1 2]);
  x = reshape(x, siz);
end
x = cat(1,x{:});
y = median(x, 1);

function y = cellmin1(x)
siz = size(x);
if siz(1)==1 && siz(2)>1
  siz([2 1]) = siz([1 2]);
  x = reshape(x, siz);
end
y = x{1};
for i=2:siz(1)
  y = min(x{i}, y);
end

function y = cellmax1(x)
siz = size(x);
if siz(1)==1 && siz(2)>1
  siz([2 1]) = siz([1 2]);
  x = reshape(x, siz);
end
y = x{1};
for i=2:siz(1)
  y = max(x{i}, y);
end

function y = celleig1(x)
% FIXME this does not work for TFR representations
siz = size(x);
if siz(1)==1 && siz(2)>1
  siz([2 1]) = siz([1 2]);
  x = reshape(x, siz);
end
x = cat(1,x{:});
% [u, s, v] = svds(real(x), 1);  % x = u * s * v'
% y = s(1,1) * v(:,1);           % retain the largest eigenvector with appropriate scaling

% this is computationally more efficient and returns a complex-valued output, following a real valued svd
[u, s] = svds(real(x*x'), 1);
y = u(:,1)'*x;

function y = cellstd1(x)
siz = size(x);
if siz(1)==1 && siz(2)>1
  siz([2 1]) = siz([1 2]);
  x = reshape(x, siz);
end
y = x{1};
n = 1;
for i=2:siz(1)
  y = y + x{i};
  n = n + 1;
end
y = std(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS to compute something over the first two dimensions of a cell-array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = cellmean2(x)
siz = size(x);
x = reshape(x, [siz(1)*siz(2) siz(3:end) 1]); % represent the first two as a single dimension
y = cellmean1(x);

function y = cellmedian2(x)
siz = size(x);
x = reshape(x, [siz(1)*siz(2) siz(3:end) 1]); % represent the first two as a single dimension
y = cellmedian1(x);

function y = cellmin2(x)
siz = size(x);
x = reshape(x, [siz(1)*siz(2) siz(3:end) 1]); % represent the first two as a single dimension
y = cellmin1(x);

function y = cellmax2(x)
siz = size(x);
x = reshape(x, [siz(1)*siz(2) siz(3:end) 1]); % represent the first two as a single dimension
y = cellmax1(x);

function y = celleig2(x)
siz = size(x);
x = reshape(x, [siz(1)*siz(2) siz(3:end) 1]); % represent the first two as a single dimension
y = celleig1(x);

function y = cellstd2(x)
siz = size(x);
x = reshape(x, [siz(1)*siz(2) siz(3:end) 1]); % represent the first two as a single dimension
y = cellstd1(x);
