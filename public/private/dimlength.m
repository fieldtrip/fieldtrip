function n = dimlength(data, seldim)

if nargin==1
  dimtok = tokenize(data.dimord, '_');
  ndim = length(dimtok);
  n = zeros(1, ndim);
  for i=1:ndim
    n(i) = dimlength(data, dimtok{i});
  end
  return
end

switch seldim
  case 'rpt'
    % hmm, this may be difficult
    if strcmp(data.dimord, 'rpt_chan_time') && isfield(data, 'trial')
      n = size(data.trial, 1);
    elseif strcmp(data.dimord, 'rpt_chan_freq') && isfield(data, 'powspctrm')
      n = size(data.powspctrm, 1);
    elseif strcmp(data.dimord, 'rpt_chan_freq') && isfield(data, 'powcovspctrm')
      n = size(data.powcovspctrm, 1);
    elseif strcmp(data.dimord, 'rpt_chan_freq') && isfield(data, 'crsspctrm')
      n = size(data.crsspctrm, 1); %in case cmbrepresentation is sparse
    elseif strcmp(data.dimord, 'rpt_chan_chan') && isfield(data, 'cov')
      n = size(data.cov, 1);
    elseif strcmp(data.dimord, 'rpt_chan_freq_time') && isfield(data, 'powspctrm')
      n = size(data.powspctrm, 1);
    elseif strcmp(data.dimord, 'rpt_chan_freq_time') && isfield(data, 'crsspctrm')
      n = size(data.crsspctrm, 1); 
    elseif strcmp(data.dimord, 'rpt_chan_chan_freq') && isfield(data, 'crsspctrm')
      n = size(data.crsspctrm, 1);
    elseif strcmp(data.dimord, 'rpt_chan_chan_freq_time') && isfield(data, 'crsspctrm')
      n = size(data.crsspctrm, 1);
    elseif strcmp(data.dimord, 'rpt_chan_freq') && isfield(data, 'fourierspctrm')
      n = size(data.fourierspctrm, 1);
    elseif strcmp(data.dimord, 'rpt_chan_freq_time') && isfield(data, 'fourierspctrm')
      n = size(data.fourierspctrm, 1);
    elseif strcmp(data.dimord, 'rpt_chan_chan_freq') && isfield(data, 'fourierspctrm')
      n = size(data.fourierspctrm, 1);
    elseif strcmp(data.dimord, 'rpt_chan_chan_freq') && isfield(data, 'powcovspctrm')
      n = size(data.powcovspctrm, 1);
    elseif strcmp(data.dimord, 'rpt_chan_chan_freq_time') && isfield(data, 'fourierspctrm')
      n = size(data.fourierspctrm, 1);
    elseif strcmp(data.dimord, 'rpt_chan_chan_freq_time') && isfield(data, 'powcovspctrm')
      n = size(data.powcovspctrm, 1);
    elseif strcmp(data.dimord, 'rpt_pos')
      %HACK to be fixed
      x = setdiff(fieldnames(data),'inside');
      for k = 1:length(x)
        dims = size(getsubfield(data,x{k}));
	if dims(2)==size(data.pos,1) && numel(dims)==2,
	  n = dims(1);
	  return
	end
      end

    elseif strcmp(data.dimord, 'rpt_pos_freq'),
      %HACK to be fixed
      x = fieldnames(data);
      for k = 1:length(x)
        dims = size(getsubfield(data,x{k}));
        if dims(2)==size(data.pos,1) && (numel(dims)==2 || dims(3)==length(data.freq)),
          n = dims(1);
	  return
	end
      end
      
    elseif strcmp(data.dimord, 'rpt_pos_time'),
      %HACK to be fixed
      x = fieldnames(data);
      for k = 1:length(x)
        dims = size(getsubfield(data,x{k}));
        if dims(2)==size(data.pos,1) && (numel(dims)==2 || dims(3)==length(data.time)),
          n = dims(1);
	  return
	end
      end
      
    else
      error('cannot determine number of repetitions for dim "%s"', seldim);
    end

  case 'rpttap'
    % hmm, this may be difficult
    if strcmp(data.dimord, 'rpttap_chan_time') && isfield(data, 'trial')
      n = size(data.trial, 1);
    elseif strcmp(data.dimord, 'rpttap_chan_freq') && isfield(data, 'powspctrm')
      n = size(data.powspctrm, 1);
    elseif strcmp(data.dimord, 'rpttap_chan_chan') && isfield(data, 'cov')
      n = size(data.cov, 1);
    elseif strcmp(data.dimord, 'rpttap_chan_freq_time') && isfield(data, 'powspctrm')
      n = size(data.powspctrm, 1);
    elseif strcmp(data.dimord, 'rpttap_chan_chan_freq') && isfield(data, 'crsspctrm')
      n = size(data.crsspctrm, 1);
    elseif strcmp(data.dimord, 'rpttap_chan_chan_freq_time') && isfield(data, 'crsspctrm')
      n = size(data.crsspctrm, 1);
    elseif strcmp(data.dimord, 'rpttap_chan_freq') && isfield(data, 'fourierspctrm')
      n = size(data.fourierspctrm, 1);
    elseif strcmp(data.dimord, 'rpttap_chan_freq_time') && isfield(data, 'fourierspctrm')
      n = size(data.fourierspctrm, 1);
    elseif strcmp(data.dimord, 'rpttap_chan_chan_freq') && isfield(data, 'fourierspctrm')
      n = size(data.fourierspctrm, 1);
    elseif strcmp(data.dimord, 'rpttap_chan_chan_freq_time') && isfield(data, 'fourierspctrm')
      n = size(data.fourierspctrm, 1);
    elseif strcmp(data.dimord, 'rpttap_pos') && isfield(data, 'crsspctrm')
      n = size(data.crsspctrm, 1);
    elseif strcmp(data.dimord, 'rpttap_pos') && isfield(data, 'fourierspctrm')
      n = size(data.fourierspctrm, 1);
    else
      error('cannot determine number of repetitions for dim "%s"', seldim);
    end

  case 'chan'
    n = length(data.label);
  case 'freq'
    n = length(data.freq);
  case 'time'
    n = length(data.time);
  case 'pos'
    n = size(data.pos,1);
  otherwise
    error('unsupported dim "%s"', seldim);
end
