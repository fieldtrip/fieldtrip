function test_ft_appendxxx(data)
% This function tests ft_appenddata handling of trialdata, particularly with
% different numbers of columns or other mismatches.  It should work with any
% data file with trialinfo and multiple channels and trials.

%% if data has multiple channels

if numel(data.label) > 1
  %% divide data into channels
  
  datasplit = cell(size(data.label));
  
  for ch = 1:numel(data.label)
    cfg           = [];
    cfg.channel   = data.label{ch};
    
    datasplit{ch} = ft_selectdata(cfg,data);
  end
  
  %% mess with some details
  
  
  
  %% remerge each channel file into one
  
  cfg = [];
  
  dataout = ft_appenddata(cfg);
  
  %% assert identity
  
  
  
end

%% if data has multiple trials

if numel(data.trial) > 1
  %% divide data into trials
  
  datasplit = cell(size(data.trial));
  
  for tr = 1:numel(data.trial)
    cfg = [];
    cfg.trial     = tr;
    
    datasplit{tr} = ft_selectdata(cfg,data);
  end
  
  %% mess with some details
  
  
  
  %% remerge each trial file into one
  
  cfg = [];
  
  dataout = ft_appenddata(cfg,datasplit);
  
  %% assert identity to input data
  
  
  
end
end