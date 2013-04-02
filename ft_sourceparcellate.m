function parcel = hcp_sourceparcellate(cfg, source)

% HCP_SOURCEPARCELLATE combines the source-reconstruction values in the
% parcels by averaging, or by concatenating.
%
% Use as
%    parcel = hcp_sourceparcellate(cfg, source)
% where the input is a 3D voxel-based or a 2D surface-based source grid
% and the output is a channel-based representation with the averages or
% concatenared representation per parcel.
%
% The configuration "cfg" is a structure that can contain the following
% fields
%   cfg.parcellation = string, fieldname that contains the desired parcellation
%   cfg.parameter    = cell-array with strings, fields that should be parcellated (can be 'all')

% Copyright (C) 2012, Robert Oostenveld

% get the defaults
cfg.parcellation = ft_getopt(cfg, 'parcellation');
cfg.parameter    = ft_getopt(cfg, 'parameter', 'all');

if ischar(cfg.parameter)
  cfg.parameter = {cfg.parameter};
end

% ensure it is a parcellation, not a segmentation
% FIXME is the indexed style needed?
source = ft_checkdata(source, 'datatype', 'parcellation', 'parcellationstyle', 'indexed');

if isempty(cfg.parcellation)
  fn = fieldnames(source);
  for i=1:numel(fn)
    if isfield(source, [fn{i} 'label'])
      warning('using %s for the parcellation', fn{i});
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
    dimord{i} = 'fixme';
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

% get the parcellation and the labels that go with it
seg      = source.(cfg.parcellation);
seglabel = source.([cfg.parcellation 'label']);
nseg     = length(seglabel);

% this will hold the output
parcel = [];
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
  univariate    = (num==1);
  bivariate     = (num==2);
  multivariate  = (num>2);
  
  if univariate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % univariate: loop over positions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear j j1 j2 % to avoid confusion
    
    for j=1:numel(seglabel)
      fprintf('parcellating %s into %s\n', fn{i}, seglabel{j});
      switch dimord{i}
        case 'pos'
          if isempty(tmp)
            tmp = nan(nseg, 1);
          end
          tmp(j) = mean(dat(seg==j), 1);
          
          
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
    
  else
    % more than two dimensions are not supported
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

try, cfg.previous = source.cfg; end
parcel.cfg = cfg;
