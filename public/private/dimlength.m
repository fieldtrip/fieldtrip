function [n, fn] = dimlength(data, seldim, fld)

fn    = fieldnames(data);
selfn = find(~cellfun('isempty', strfind(fn, 'dimord')));
fn    = fn(selfn);

for k = 1:numel(fn)
  fndimord{k} = data.(fn{k});
end

if nargin==1
  for k = 1:numel(fn)
    dimtok = tokenize(fndimord{k}, '_');
    ndim = length(dimtok);
    n{k,1} = zeros(1, ndim);
    for i=1:ndim
      n{k}(i) = dimlength(data, dimtok{i}, fn{k});
    end
  end
  return
end

switch seldim
  case 'rpt'
    % hmm, this may be difficult
    if numel(fld)>6 && isfield(data, fld(1:end-6)),
      dimtok = tokenize(data.(fld), '_'); 
      tmp    = data.(fld(1:end-6));
      if iscell(tmp)
        if isfield(data, 'inside'),
          ix = data.inside(1);
        else
          ix = 1;
        end
        tmp    = tmp{ix};
        dimtok = dimtok(2:end);
      end
      ix = find(~cellfun('isempty', strfind(dimtok, seldim)));
      n  = size(tmp, ix);
    elseif strcmp(data.(fld), 'rpt_chan_time') && isfield(data, 'trial')
      n = size(data.trial, 1);
    elseif strcmp(data.(fld), 'rpt_chan_freq') && isfield(data, 'powspctrm')
      n = size(data.powspctrm, 1);
    elseif strcmp(data.(fld), 'rpt_chan_freq') && isfield(data, 'powcovspctrm')
      n = size(data.powcovspctrm, 1);
    elseif strcmp(data.(fld), 'rpt_chan_freq') && isfield(data, 'crsspctrm')
      n = size(data.crsspctrm, 1); %in case cmbrepresentation is sparse
    elseif strcmp(data.(fld), 'rpt_chan_chan') && isfield(data, 'cov')
      n = size(data.cov, 1);
    elseif strcmp(data.(fld), 'rpt_chan_freq_time') && isfield(data, 'powspctrm')
      n = size(data.powspctrm, 1);
    elseif strcmp(data.(fld), 'rpt_chan_freq_time') && isfield(data, 'crsspctrm')
      n = size(data.crsspctrm, 1); 
    elseif strcmp(data.(fld), 'rpt_chan_chan_freq') && isfield(data, 'crsspctrm')
      n = size(data.crsspctrm, 1);
    elseif strcmp(data.(fld), 'rpt_chan_chan_freq_time') && isfield(data, 'crsspctrm')
      n = size(data.crsspctrm, 1);
    elseif strcmp(data.(fld), 'rpt_chan_freq') && isfield(data, 'fourierspctrm')
      n = size(data.fourierspctrm, 1);
    elseif strcmp(data.(fld), 'rpt_chan_freq_time') && isfield(data, 'fourierspctrm')
      n = size(data.fourierspctrm, 1);
    elseif strcmp(data.(fld), 'rpt_chan_chan_freq') && isfield(data, 'fourierspctrm')
      n = size(data.fourierspctrm, 1);
    elseif strcmp(data.(fld), 'rpt_chan_chan_freq') && isfield(data, 'powcovspctrm')
      n = size(data.powcovspctrm, 1);
    elseif strcmp(data.(fld), 'rpt_chan_chan_freq_time') && isfield(data, 'fourierspctrm')
      n = size(data.fourierspctrm, 1);
    elseif strcmp(data.(fld), 'rpt_chan_chan_freq_time') && isfield(data, 'powcovspctrm')
      n = size(data.powcovspctrm, 1);
    elseif strcmp(data.(fld), 'rpt_pos')
      %HACK to be fixed
      x = setdiff(fld(data),'inside');
      for k = 1:length(x)
        dims = size(getsubfield(data,x{k}));
	if dims(2)==size(data.pos,1) && numel(dims)==2,
	  n = dims(1);
	  return
	end
      end

    elseif strcmp(data.(fld), 'rpt_pos_freq'),
      %HACK to be fixed
      x = fld(data);
      for k = 1:length(x)
        dims = size(getsubfield(data,x{k}));
        if dims(2)==size(data.pos,1) && (numel(dims)==2 || dims(3)==length(data.freq)),
          n = dims(1);
	  return
	end
      end
      
    elseif strcmp(data.(fld), 'rpt_pos_time'),
      %HACK to be fixed
      x = fld(data);
      for k = 1:length(x)
        dims = size(getsubfield(data,x{k}));
        if dims(2)==size(data.pos,1) && (numel(dims)==2 || dims(3)==length(data.time)),
          n = dims(1);
	  return
	end
      end
      
    else
      %error('cannot determine number of repetitions for dim "%s"', seldim);
      n = nan;
    end

  case 'rpttap'
    % hmm, this may be difficult
    if numel(fld)>6 && isfield(data, fld(1:end-6)),
      dimtok = tokenize(data.(fld), '_'); 
      tmp    = data.(fld(1:end-6));
      if iscell(tmp)
        if isfield(data, 'inside'),
          ix = data.inside(1);
        else
          ix = 1;
        end
        tmp    = tmp{ix};
        dimtok = dimtok(2:end);
      end
      ix = find(~cellfun('isempty', strfind(dimtok, seldim)));
      n  = size(tmp, ix);
    elseif strcmp(data.(fld), 'rpttap_chan_time') && isfield(data, 'trial')
      n = size(data.trial, 1);
    elseif strcmp(data.(fld), 'rpttap_chan_freq') && isfield(data, 'powspctrm')
      n = size(data.powspctrm, 1);
    elseif strcmp(data.(fld), 'rpttap_chan_chan') && isfield(data, 'cov')
      n = size(data.cov, 1);
    elseif strcmp(data.(fld), 'rpttap_chan_freq_time') && isfield(data, 'powspctrm')
      n = size(data.powspctrm, 1);
    elseif strcmp(data.(fld), 'rpttap_chan_chan_freq') && isfield(data, 'crsspctrm')
      n = size(data.crsspctrm, 1);
    elseif strcmp(data.(fld), 'rpttap_chan_chan_freq_time') && isfield(data, 'crsspctrm')
      n = size(data.crsspctrm, 1);
    elseif strcmp(data.(fld), 'rpttap_chan_freq') && isfield(data, 'fourierspctrm')
      n = size(data.fourierspctrm, 1);
    elseif strcmp(data.(fld), 'rpttap_chan_freq_time') && isfield(data, 'fourierspctrm')
      n = size(data.fourierspctrm, 1);
    elseif strcmp(data.(fld), 'rpttap_chan_chan_freq') && isfield(data, 'fourierspctrm')
      n = size(data.fourierspctrm, 1);
    elseif strcmp(data.(fld), 'rpttap_chan_chan_freq_time') && isfield(data, 'fourierspctrm')
      n = size(data.fourierspctrm, 1);
    elseif strcmp(data.(fld), 'rpttap_pos') && isfield(data, 'crsspctrm')
      n = size(data.crsspctrm, 1);
    elseif strcmp(data.(fld), 'rpttap_pos') && isfield(data, 'fourierspctrm')
      n = size(data.fourierspctrm, 1);
    else
      %error('cannot determine number of repetitions for dim "%s"', seldim);
      n = nan;
    end

  case 'chan'
    if ~isfield(data, 'inside'), 
      try, 
        n = length(data.label);
      catch
        n = size(data.labelcmb, 1);
      end
    else
      n = nan; %FIXME discuss appending label to source-like data
    end
  case 'freq'
    n = length(data.freq);
  case 'time'
    n = length(data.time);
  case {'pos' '{pos}'}
    n = size(data.pos,1);
  case {'ori'}
    if isfield(data, 'ori'), 
      n = size(data.ori,1);
    else
      n = 1;
    end
  otherwise
    error('unsupported dim "%s"', seldim);
end
