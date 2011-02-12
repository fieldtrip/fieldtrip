function [source] = fixsource(source)

% FIXSOURCE converts old style source structures into new style source structures
%
% Use as
%   output = fixsource(input)
% where input is a structure representing source data
%
% Typically, old style source structures contain source.avg.XXX or source.trial.XXX. 
% Furthermore. old style source structrures do not contain a dimord field.
%
% The new style source structure contains:
%   source.pos       Nx3 list with source positions
%   source.dim       optional, if the list of positions describes a 3D volume
%   source.XXX       the old style subfields in avg/trial
%   source.XXXdimord string how to interpret the respective XXX field:
%
% For example 
%   source.pow             = zeros(Npos,Ntime)
%   source.powdimord       = 'pos_time'
%
%   source.mom             = cell(1,Npos)
%   source.mom{1}          = zeros(Nori,Nrpttap)
%   source.momdimord       = '{pos}_ori_rpttap'
%
%   source.leadfield       = cell(1,Npos)
%   source.leadfield{1}    = zeros(Nchan,Nori)
%   source.leadfielddimord = '{pos}_chan_ori'
%
% See also FT_CHECKDATA, FIXVOLUME

npos = size(source.pos, 1);

if ~isfield(source, 'inside')
  source.inside = 1:npos;
  source.outside = [];
end

% old data structures may use latency/frequency instead of time/freq
if isfield(data, 'frequency'),
  data.freq = data.frequency;
  data      = rmfield(data, 'frequency');
end
if isfield(data, 'latency'),
  data.time = data.latency;
  data      = rmfield(data, 'latency');
end

% copy the fields from the average into the main structure
if isfield(source, 'avg')
  avg    = source.avg;
  source = rmfield(source,  'avg');
  
  fn = fieldnames(avg);
  for i=1:length(fn)
    source.(fn{i}) = avg.(fn{i});
  end
end

% convering the trials has precedence over the average
if isfield(source, 'trial')
  trial  = source.trial;
  source = rmfield(source,  'trial');
  nrpt   = length(trial);
  for i=1:nrpt
    % FIXME
    keyboard
  end
end

hastime = isfield(source, 'time');
hasfreq = isfield(source, 'freq');

if hastime
  ntime = length(source.time);
else
  ntime = 1;
end

if hasfreq
  nfreq = length(source.freq);
else
  nfreq = 1;
end


% FIXME the code bellow should be repatterned
% first get the dimord string
% then determine the data dimensions
% then determine the dimensions according to the dimord
% if consistent, continue
% if inconsistent but can be reshaped, reshape and give warning
% if inconsistent and cannot can be reshaped, give error


% check the dimord and dimensions for each of the fields
fn = fieldnames(source);
for i=1:length(fn)
  element = source.(fn{i});
  
  if ischar(element) || isstruct(element)
    % the dimord only applies to cell-arrays or numeric arrays
    continue
  end
  
  if any(strcmp(fn{i}, {'inside', 'outside'}))
    % these don't need a dimord
    continue
  end

  % determine the dimensions of this data element
  if iscell(element)
    dim = numel(element);
    dim = [dim size(element{source.inside(1)})];
  else
    dim = size(element);
  end
  
  if isfield(source, [fn{i} 'dimord'])
    % check whether the specific dimord applies
    dimord = source.([fn{i} 'dimord']);
    
    dimtok = tokenize(dimord, '_');
    descr  = [];
    for j=1:length(dimtok)
      if strcmp(dimtok{j}, 'pos')
        descr = [descr size(source.pos,1)];
      elseif strcmp(dimtok{j}, '{pos}')
        descr = [descr size(source.pos,1)];
      else
        descr = [descr numel(source.(dimtok{j}))];
      end
    end
    
    if prod(dim)~=prod(descr)
      error('the dimensions of the source data are not consistent with the dimord (%s)', fn{i})
    end
    
    if iscell(element)
      element = reshape(element, [descr(1) 1]);
      for j=1:numel(element)
        element{j} = reshape(element{j}, [descr(2:end) 1]);
      end
    else
      element = reshape(element, [descr 1]);
    end
    
    % put the reshaped data element back into the structure
    source.(fn{i}) = element;
    
  elseif isfield(source, 'dimord')
    % check whether the general dimord applies to this element
    dimord = source.dimord;
    
    dimtok = tokenize(dimord, '_');
    descr  = [];
    for j=1:length(dimtok)
      if strcmp(dimtok{j}, 'pos')
        descr = [descr size(source.pos,1)];
      elseif strcmp(dimtok{j}, '{pos}')
        descr = [descr size(source.pos,1)];
      else
        descr = [descr numel(source.(dimtok{j}))];
      end
    end
    
    if prod(dim)==prod(descr)
      % the dimensions of the source data are consistent with the dimord
      if iscell(element)
        element = reshape(element, [descr(1) 1]);
        for j=1:numel(element)
          element{j} = reshape(element{j}, [descr(2:end) 1]);
        end
      else
        element = reshape(element, [descr 1]);
      end
      
      % create a dimord for this data element
      source.([fn{i} 'dimord']) = dimord;
      % put the reshaped data element back into the structure
      source.(fn{i}) = element;
    end
    
  else
    % try to create a dimord on the fly
    if     ~hasfreq && ~hastime && prod(dim)==npos
      if iscell(element)
        dimord = '{pos}';
      else
        dimord = 'pos';
      end
      descr  = npos;
    elseif  hasfreq && ~hastime && prod(dim)==prod([npos nfreq])
      if iscell(element)
        dimord = '{pos}_freq';
      else
        dimord = 'pos_freq';
      end
      descr  = [npos nfreq];
    elseif  hasfreq &&  hastime && prod(dim)==prod([npos nfreq ntime])
      if iscell(element)
        dimord = '{pos}_freq_time';
      else
        dimord = 'pos_freq_time';
      end
      descr  = [npos nfreq ntime];
    elseif ~hasfreq &&  hastime && prod(dim)==prod([npos ntime])
      if iscell(element)
        dimord = '{pos}_time';
      else
        dimord = 'pos_time';
      end
      descr  = [npos ntime];
    elseif prod(dim)==npos
      % although the source data has time and frequency, this element might be pos only
      if iscell(element)
        dimord = '{pos}';
      else
        dimord = 'pos';
      end
      descr  = npos;
    else
      dimord = '';
      descr  = [];
    end
    
    if ~isempty(dimord)
      % the dimensions of the source data are consistent with the dimord
      if iscell(element)
        element = reshape(element, [descr(1) 1]);
        for j=1:numel(element)
          element{j} = reshape(element{j}, [descr(2:end) 1]);
        end
      else
        element = reshape(element, [descr 1]);
      end
      
      % create a dimord for this data element
      source.([fn{i} 'dimord']) = dimord;
      % put the reshaped data element back into the structure
      source.(fn{i}) = element;
    end
    
  end % if isfield dimord
end % for each of the fields


if isfield(source, 'dimord')
  % remove the general dimord field, there should be one for each data element
  source = rmfield(source, 'dimord');
end


%if isvolume
%  % ensure consistent dimensions of the volumetric data
%  % reshape each of the volumes that is found into a 3D array
%  param = parameterselection('all', data);
%  dim   = data.dim;
%  for i=1:length(param)
%    tmp  = getsubfield(data, param{i});
%    tmp  = reshape(tmp, dim);
%    data = setsubfield(data, param{i}, tmp);
%  end
%end


% if issource || isvolume,
%   % the following section is to make a dimord-consistent representation of
%   % volume and source data, taking trials, time and frequency into account
%   if isequal(hasdimord, 'yes') && (~isfield(data, 'dimord') || ~strcmp(data.dimord, sourcedimord))
%
%     % determine the size of the data
%     if isfield(data, 'dimord'),
%       dimtok = tokenize(data.dimord, '_');
%       if ~isempty(strmatch('time', dimtok)), Ntime = length(data.time); else Ntime = 1; end
%       if ~isempty(strmatch('freq', dimtok)), Nfreq = length(data.freq); else Nfreq = 1; end
%     else
%       Nfreq = 1;
%       Ntime = 1;
%     end
%
%     % convert old style source representation into new style
%     if isfield(data, 'avg') && isfield(data.avg, 'mom') && isfield(data, 'freq') && strcmp(sourcedimord, 'rpt_pos'),
%       % frequency domain source representation convert to single trial power
%       Npos   = size(data.pos,1);
%       Nrpt   = length(data.cumtapcnt);
%       tmpmom = zeros(Npos, size(data.avg.mom{data.inside(1)},2));
%       tmpmom(data.inside,:) = cat(1,data.avg.mom{data.inside});
%       tmppow = zeros(Npos, Nrpt);
%       tapcnt = [0;cumsum(data.cumtapcnt)];
%       for k = 1:Nrpt
%         Ntap = tapcnt(k+1)-tapcnt(k);
%         tmppow(data.inside,k) = sum(abs(tmpmom(data.inside,(tapcnt(k)+1):tapcnt(k+1))).^2,2)./Ntap;
%       end
%       data.pow = tmppow';
%       data     = rmfield(data, 'avg');
%       if strcmp(inside, 'logical'),
%         data     = fixinside(data, 'logical');
%         data.inside = repmat(data.inside(:)',[Nrpt 1]);
%       end
%     elseif isfield(data, 'avg') && isfield(data.avg, 'mom') && isfield(data, 'freq') && strcmp(sourcedimord, 'rpttap_pos'),
%       % frequency domain source representation convert to single taper fourier coefficients
%       Npos   = size(data.pos,1);
%       Nrpt   = sum(data.cumtapcnt);
%       data.fourierspctrm = complex(zeros(Nrpt, Npos), zeros(Nrpt, Npos));
%       data.fourierspctrm(:, data.inside) = transpose(cat(1, data.avg.mom{data.inside}));
%       data   = rmfield(data, 'avg');
%     elseif isfield(data, 'avg') && isfield(data.avg, 'mom') && isfield(data, 'time') && strcmp(sourcedimord, 'pos_time'),
%       Npos   = size(data.pos,1);
%       Nrpt   = 1;
%       tmpmom = zeros(Npos, size(data.avg.mom{data.inside(1)},2));
%       tmpmom(data.inside,:) = cat(1,data.avg.mom{data.inside});
%       data.mom = tmpmom;
%       if isfield(data.avg, 'noise'),
%         tmpnoise = data.avg.noise(:);
%         data.noise = tmpnoise(:,ones(1,size(tmpmom,2)));
%       end
%       data = rmfield(data, 'avg');
%       Ntime = length(data.time);
%     elseif isfield(data, 'trial') && isfield(data.trial(1), 'mom') && isfield(data, 'time') && strcmp(sourcedimord, 'rpt_pos_time'),
%       Npos   = size(data.pos,1);
%       Nrpt   = length(data.trial);
%       Ntime  = length(data.time);
%       tmpmom = zeros(Nrpt, Npos, Ntime);
%       for k = 1:Nrpt
%         tmpmom(k,data.inside,:) = cat(1,data.trial(k).mom{data.inside});
%       end
%       data     = rmfield(data, 'trial');
%       data.mom = tmpmom;
%     elseif isfield(data, 'trial') && isstruct(data.trial)
%       Nrpt = length(data.trial);
%     else
%       Nrpt = 1;
%     end
%
%     % start with an initial specification of the dimord and dim
%     if (~isfield(data, 'dim') || ~isfield(data, 'dimord'))
%       if issource
%         % at least it should have a Nx3 pos
%         data.dim    = size(data.pos, 1);
%         data.dimord = 'pos';
%       elseif isvolume
%         % at least it should have a 1x3 dim
%         data.dim    = data.dim;
%         data.dimord = 'dim1_dim2_dim3';
%       end
%     end
%
%     % add the additional dimensions
%     if Nfreq>1
%       data.dimord = [data.dimord '_freq'];
%       data.dim    = [data.dim     Nfreq];
%     end
%     if Ntime>1
%       data.dimord = [data.dimord '_time'];
%       data.dim    = [data.dim     Ntime];
%     end
%     if Nrpt>1 && strcmp(sourcedimord, 'rpt_pos'),
%       data.dimord = ['rpt_' data.dimord];
%       data.dim    = [Nrpt   data.dim ];
%     elseif Nrpt>1 && strcmp(sourcedimord, 'rpttap_pos'),
%       data.dimord = ['rpttap_' data.dimord];
%       data.dim    = [Nrpt   data.dim ];
%     end
%
%     % the nested trial structure is not compatible with dimord
%     if isfield(data, 'trial') && isstruct(data.trial)
%       param = fieldnames(data.trial);
%       for i=1:length(param)
%         if isa(data.trial(1).(param{i}), 'cell')
%           concat = cell(data.dim(1), prod(data.dim(2:end)));
%         else
%           concat = zeros(data.dim(1), prod(data.dim(2:end)));
%         end
%         for j=1:length(data.trial)
%           tmp = data.trial(j).(param{i});
%           concat(j,:) = tmp(:);
%         end % for each trial
%         data.trial = rmfield(data.trial, param{i});
%         data.(param{i}) = reshape(concat, data.dim);
%       end % for each param
%       data = rmfield(data, 'trial');
%     end
%   end
%
%   % ensure consistent dimensions of the source reconstructed data
%   % reshape each of the source reconstructed parameters
%   if issource && isfield(data, 'dim') && prod(data.dim)==size(data.pos,1)
%     % volumetric represented source structure
%     dim = [prod(data.dim) 1];
%   elseif issource && any(~cellfun('isempty',strfind(fieldnames(data), 'dimord')))
%     % sparsely represented source structure new style
%     dim = [size(data.pos,1) 1];
%   elseif isfield(data, 'dim')
%     dim = [data.dim 1];
%   elseif isfield(data, 'dimord')
%     %HACK
%     dimtok = tokenize(data.dimord, '_');
%     for i=1:length(dimtok)
%       if strcmp(dimtok(i), 'pos')
%         dim(1,i) = size(getsubfield(data,dimtok{i}),1);
%       elseif strcmp(dimtok(i), 'rpt')
%         dim(1,i) = nan;
%       else
%         dim(1,i) = length(getsubfield(data,dimtok{i}));
%       end
%     end
%     i = find(isnan(dim));
%     if ~isempty(i)
%       n = fieldnames(data);
%       for ii=1:length(n)
%         numels(1,ii) = numel(getfield(data,n{ii}));
%       end
%       nrpt = numels./prod(dim(setdiff(1:length(dim),i)));
%       nrpt = nrpt(nrpt==round(nrpt));
%       dim(i) = max(nrpt);
%     end
%     if numel(dim)==1, dim(1,2) = 1; end;
%   end
%
% end
