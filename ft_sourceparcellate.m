function parcel = ft_sourceparcellate(cfg, source, parcellation)

% FT_SOURCEPARCELLATE combines the source-reconstruction values in the
% parcels by averaging or by concatenating.
%
% Use as
%    parcel = ft_sourceparcellate(cfg, source)
% where the input is a 3D voxel-based or a 2D surface-based source grid
% and the output is a channel-based representation with the averages or
% concatenared representation per parcel.
%
% The input data structure is for example obtained from FT_SOURCEANALYSIS
% or FT_COMPUTE_LEADFIELD and is described in detail in FT_DATATYPE_SOURCE.
% The output data structure is described in detail in
% FT_DATATYPE_PARCELLATION.
%
% The configuration "cfg" is a structure that can contain the following
% fields
%   cfg.method       = string
%   cfg.parcellation = string, fieldname that contains the desired parcellation
%   cfg.parameter    = cell-array with strings, fields that should be parcellated (default = 'all')
%
% See also FT_SOURCEANALYSIS, FT_DATATYPE_PARCELLATION, FT_DATATYPE_SEGMENTATION

% Copyright (C) 2012-2013, Robert Oostenveld
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

revision = '$Id$';

ft_defaults
ft_preamble help
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug
ft_preamble loadvar source

% get the defaults
cfg.parcellation = ft_getopt(cfg, 'parcellation');
cfg.parameter    = ft_getopt(cfg, 'parameter', 'all');

if ischar(cfg.parameter)
  cfg.parameter = {cfg.parameter};
end

if nargin<3
  % the parcellation is represented in the source structure
  parcellation = source;
else
  % the parcellation is specified as separate structure
end

% ensure it is a parcellation, not a segmentation
parcellation = ft_checkdata(parcellation, 'datatype', 'parcellation', 'parcellationstyle', 'indexed'); % FIXME is the indexed style needed?
source       = ft_checkdata(source, 'datatype', 'source', 'inside', 'logical', 'sourcerepresentation', 'new');

if isempty(cfg.parcellation)
  fn = fieldnames(parcellation);
  for i=1:numel(fn)
    if isfield(parcellation, [fn{i} 'label'])
      warning('using "%s" for the parcellation', fn{i});
      cfg.parcellation = fn{i};
      break
    end
  end
end

if isempty(cfg.parcellation)
  error('you should specify the field containing the parcellation');
end

if isfield(source, 'dimord')
  % determine the size of fields that are consistent with the general dimord
  tok = tokenize(source.dimord, '_');
  siz = nan(size(tok));
  for i=1:length(tok)
    switch tok{i}
      case 'pos'
        siz(i) = size(source.pos,1);
      case 'time'
        siz(i) = length(source.time);
      case 'freq'
        siz(i) = length(source.freq);
      otherwise
        error('cannot determine the dimensions for "%s"', tok{i});
    end % switch
  end
  if numel(siz)==1
    % the size function always returns 2 or more elements
    siz(2) = 1;
  end
else
  % this will cause a failure further down in the code
  siz = nan;
end

fn     = fieldnames(source);
sel    = false(size(fn));
dimord = cell(size(fn));
for i=1:numel(fn)
  tmp = source.(fn{i});
  if isfield(source, [fn{i} 'dimord'])
    sel(i)    = true;
    dimord{i} = source.([fn{i} 'dimord']); % a specific dimord
  elseif iscell(tmp) && numel(tmp)==size(source.pos,1)
    sel(i)    = true;
    dimord{i} = '{pos}';
  elseif ~iscell(tmp) && isequal(size(tmp), siz)
    sel(i)    = true;
    dimord{i} = source.dimord; % the general dimord
  end
end

% these two will now contain the fields and corresponding dimord to work on
fn     = fn(sel);
dimord = dimord(sel);

if any(strcmp(cfg.parameter, 'all'))
  cfg.parameter = fn;
else
  [inside, i1, i2] = intersect(cfg.parameter, fn);
  [outside       ] = setdiff(cfg.parameter, fn);
  if ~isempty(outside)
    warning('\nparameter "%s" cannot be parcellated', outside{:});
  end
  cfg.parameter = fn(i2);
  fn     = fn(i2);
  dimord = dimord(i2);
end

if numel(fn)==0
  error('there are no source parameters that can be parcellated');
end

% get the parcellation and the labels that go with it
seg      = parcellation.(cfg.parcellation);
seglabel = parcellation.([cfg.parcellation 'label']);
nseg     = length(seglabel);

if isfield(source, 'inside')
  % determine the conjunction of the parcellation and the inside source points
  % points that are
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

% this will hold the output
parcel       = [];
parcel.label = seglabel;
if isfield(source, 'time')
  parcel.time = source.time;
  ntime = length(source.time);
else
  ntime = 1; % for convenience
end
if isfield(source, 'freq')
  parcel.freq = source.freq;
  nfreq = length(source.freq);
else
  nfreq = 1; % for convenience
end

for i=1:numel(fn)
  if strcmp(fn{i}, cfg.parcellation)
    % don't parcellate the parcellation itself
    continue
  end
  
  dat = source.(fn{i});
  tmp = [];
  
  % determine how many "position" dimensions there are
  num = numel(strmatch('pos', tokenize(dimord{i}, '_')));
  if num==0
    % the positions might be represented in a cell-array, i.e. {pos}_time
    num = numel(strmatch('{pos}', tokenize(dimord{i}, '_')));
  end
  univariate    = (num==1);
  bivariate     = (num==2);
  multivariate  = (num>2);
  
  if univariate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % univariate: loop over positions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%]
    clear j j1 j2 % to avoid confusion
    
    for j=1:numel(seglabel)
      fprintf('parcellating %s into %s\n', fn{i}, seglabel{j});
      switch dimord{i}
        case 'pos'
          if isempty(tmp)
            tmp = nan(nseg, 1);
          end
          tmp(j) = mean(dat(seg==j), 1);
          
        case {'{pos}' '{pos}_ori_ori', '{pos}_ori_time'}
          if isempty(tmp)
            tmp = cell(nseg, 1);
          end
          tmp{j} = cellmean(dat(seg==j));
          
        case 'pos_time'
          if isempty(tmp)
            tmp = nan(nseg, ntime);
          end
          tmp(j,:) = mean(dat(seg==j,:), 1);
          
        case 'pos_freq'
          if isempty(tmp)
            tmp = nan(nseg, nfreq);
          end
          tmp(j,:) = mean(dat(seg==j,:), 1);
          
        case 'pos_freq_time'
          if isempty(tmp)
            tmp = nan(nseg, nfreq, ntime);
          end
          tmp(j,:,:) = mean(dat(seg==j,:,:), 1);
          
          % case '{pos}'
          % case '{pos}_time'
          % case '{pos}'
          % case '{pos}_freq'
          % etc.
          
        otherwise
          error('dimord "%s" cannot be parcellated', dimord{i});
      end % switch
    end % for each of the parcels
    
    
  elseif bivariate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % bivariate: loop over pairs of positions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear j j1 j2 % to avoid confusion
    
    for j1=1:numel(seglabel)
      for j2=1:numel(seglabel)
        fprintf('parcellating %s into %s\n', fn{i}, seglabel{j});
        switch dimord{i}
          case 'pos_pos'
            if isempty(tmp)
              tmp = nan(nseg, nseg);
            end
            tmp(j1,j2) = mean(mean(dat(seg==j1,seg==j2), 1), 2);
            
          case 'pos_pos_time'
            if isempty(tmp)
              tmp = nan(nseg, nseg, ntime);
            end
            tmp(j1,j2,:) = mean(mean(dat(seg==j1,seg==j2,:), 1), 2);
            
          case 'pos_pos_freq'
            if isempty(tmp)
              tmp = nan(nseg, nseg, nfreq);
            end
            tmp(j1,j2,:) = mean(mean(dat(seg==j1,seg==j2,:), 1), 2);
            
          case 'pos_pos_freq_time'
            if isempty(tmp)
              tmp = nan(nseg, nseg, nfreq, ntime);
            end
            tmp(j1,j2,:,:) = mean(mean(dat(seg==j1,seg==j2,:,:), 1), 2);
            
            % case '{pos_pos}_time'
            % case '{pos_pos}_freq'
            % etc.
            
          otherwise
            error('dimord "%s" cannot be parcellated', dimord{i});
        end % switch
      end % for each of the parcels
    end % for each of the parcels
    
  elseif multivariate
    % more than two dimensions are not supported
    error('dimord "%s" cannot be parcellated', dimord{i});
    
  else
    error('dimord "%s" cannot be parcellated', dimord{i});
  end
  
  % update the dimord, use chan rather than pos
  % this makes it look just like timelock or freq data
  tok = tokenize(dimord{i}, '_');
  tok(strcmp(tok,  'pos' )) = { 'chan' }; % replace pos by chan
  tok(strcmp(tok, '{pos}')) = {'{chan}'}; % replace pos by chan
  tmpdimord = sprintf('%s_', tok{:});
  tmpdimord = tmpdimord(1:end-1); % exclude the last _
  
  % store the results in the output structure
  parcel.(fn{i})            = tmp;
  parcel.([fn{i} 'dimord']) = tmpdimord;
  
  clear dat tmp tmpdimord
end % for each of the fields that should be parcellated

ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble previous source parcellation
ft_postamble history parcel
ft_postamble savevar parcel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = cellmean(x)
siz = size(x);
if sum(siz==1)>1
  error('invalid size of cell array');
end
y = x{1};
for i=2:prod(siz)
  y = y + x{i};
end
y = y/prod(siz);

