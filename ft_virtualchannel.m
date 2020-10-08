function [parcel] = ft_virtualchannel(cfg, data, source, parcellation)

% FT_VIRTUALCHANNEL combines data, that is represented at the level of 
% individual dipole locations, across dipoles that belong to the same 
% anatomically of functionally parcel, for instance by averaging.
%
% Use as
%    output = ft_sourceparcellate(cfg, source)
% or 
%    output = ft_sourceparcellate(cfg, source, parcellation)
% or 
%    output = ft_sourceparcellate(cfg, source, data)
% 
% where the input source is a 2-D surface-based or 3-D voxel-based source
% grid that was for example obtained from FT_SOURCEANALYSIS or 
% FT_PREPARE_LEADFIELD. If the function is used with only a single input
% data argument, or with an additional data argument that is defined at the 
% channel level, the variable source should also contain information about
% parcellation scheme. Alternatively, the parcellation to be used can be a
% separate input argument. The input parcellation is described in detail in
% FT_DATATYPE_PARCELLATION (2-D) or FT_DATATYPE_SEGMENTATION (3-D) and can 
% be obtained from FT_READ_ATLAS or from a custom parcellation/segmentation
% for your individual subject. If the function is used with an additional
% channel level data argument, the input source should contain the spatial
% filters that map from sensor to source level data.
% The output is a channel-based data representation with the combined (e.g.
% averaged) representation of the specified parameters per parcel. In the
% presence of a data input argument, the parameters should refer to this
% argument, otherwise the parameters will be selected from the source
% input.
%
% The configuration "cfg" is a structure that can contain the following
% fields
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

% Copyright (C) 2012-2019, Robert Oostenveld
% Copyright (C) 2020, Jan-Mathijs Schoffelen and Robert Oostenveld
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
  siz = getdimsiz(source, fn{i});
  siz(contains(tokenize(dimord{i},'_'),'pos')) = nseg;
  
  if strncmp('{pos_pos}', dimord{i}, 9)
    fprintf('creating %d*%d parcel combinations for parameter %s by taking the %s\n', numel(seglabel), numel(seglabel), fn{i}, cfg.method);
    tmp = zeros(siz);
    ft_progress('init', cfg.feedback, 'computing parcellation');
    k = 0;
    K = numel(seglabel)^2;
    for j1=1:numel(seglabel)
      for j2=1:numel(seglabel)
        k = k + 1;
        ft_progress(k/K, 'computing parcellation for %s combined with %s', seglabel{j1}, seglabel{j2});
        switch cfg.method
          case 'mean'
            tmp(j1,j2,:,:) = cellmean2(dat(seg==j1,seg==j2,:));
          case 'median'
            tmp(j1,j2,:,:) = cellmedian2(dat(seg==j1,seg==j2,:));
          case 'min'
            tmp(j1,j2,:,:) = cellmin2(dat(seg==j1,seg==j2,:));
          case 'max'
            tmp(j1,j2,:,:) = cellmax2(dat(seg==j1,seg==j2,:));
          case 'eig'
            tmp(j1,j2,:,:) = celleig2(dat(seg==j1,seg==j2,:));
          case 'std'
            tmp(j1,j2,:,:) = cellstd2(dat(seg==j1,seg==j2,:));
          otherwise
            ft_error('method %s not implemented for %s', cfg.method, dimord{i});
        end % switch
      end % for j2
    end % for j1
    ft_progress('close');
    
  elseif strncmp('{pos}', dimord{i}, 5)
    fprintf('creating %d parcels for parameter %s by taking the %s\n', numel(seglabel), fn{i}, cfg.method);
    tmp = zeros(siz);
    ft_progress('init', cfg.feedback, 'computing parcellation');
    for j=1:numel(seglabel)
      ft_progress(j/numel(seglabel), 'computing parcellation for %s', seglabel{j});
      switch cfg.method
        case 'mean'
          tmp(j,:,:) = cellmean1(dat(seg==j));
        case 'median'
          tmp(j,:,:) = cellmedian1(dat(seg==j));
        case 'min'
          tmp(j,:,:) = cellmin1(dat(seg==j));
        case 'max'
          tmp(j,:,:) = cellmax1(dat(seg==j));
        case 'eig'
          tmp(j,:,:) = celleig1(dat(seg==j));
        case 'std'
          tmp(j,:,:) = cellstd1(dat(seg==j));
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

ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   source parcellation
ft_postamble provenance parcel
ft_postamble history    parcel
ft_postamble savevar    parcel
