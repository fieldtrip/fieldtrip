function test_ft_appenddata(data,mismatch)
% This function tests ft_appenddata handling of trialdata, particularly with
% different numbers of columns or other mismatches.
%
%   INPUTS:
%     data      = any FT data structure with trialinfo and multiple channels
%         and/or trials - these will be split and re-merged by ft_appenddata
%     mismatch  = boolean flag for whether you want the function to artificially
%         introduce a mismatch in trialinfo (TRUE, default) and test whether it
%         is successfully identified and removed from the output, or make only
%         inconsequential changes to trialinfo (FALSE) and verify that the
%         append function runs successfully, retaining the correct trialinfo
%
% written 3/16/17 by Teresa E. Madsen

%% check inputs

if ~isfield(data,'trialinfo')
  error('data.trialinfo does not exist, so there is nothing to test!')
end

if nargin < 2 || isempty(mismatch)
  mismatch = true;
else
  mismatch = istrue(mismatch);
end

%% if data has multiple channels, they will be split & re-merged within trials
% so trialinfo must be identical

if numel(data.label) > 1
  %% divide data into channels
  
  cfg           = [];
  cfg.channel   = data.label{1};
  data1 = ft_selectdata(cfg,data);
  
  cfg.channel   = data.label{2:end};
  data2 = ft_selectdata(cfg,data);
  
  %% mess with some details
  
  if mismatch
    data1.trialinfo(1,1) = data1.trialinfo(1,1) + 1;
  end
  
  %% remerge each channel file into one
  
  dataout = ft_appenddata([],data1,data2);
  
  %% assert identity
  
  if mismatch
    assert(~isfield(dataout,'trialinfo'),'mismatched trialinfo was not removed')
    disp('mismatched trialinfo was successfully removed upon merge by channels')
  else
    assert(isequaln(data.trialinfo,dataout.trialinfo),'something changed from input to output data')
    % data.hdr not included, data.cfg should change
    disp('trialinfo passed through split & merge by channels unchanged')
  end
end

%% if data has multiple trials, they will be split & re-merged across trials
% so trialinfo just has to have the same # of columns

if numel(data.trial) > 1
  %% separate data by trials
  
  cfg = [];
  cfg.trials = 1;
  data1 = ft_selectdata(cfg,data);
  
  cfg.trials = 2:numel(data.trial);
  data2 = ft_selectdata(cfg,data);
  
  %% mess with some details
  
  if mismatch   % adding a column prevents concatenation across trials
    data1.trialinfo = [data1.trialinfo zeros(1,1)];
  else  % this should be okay, as these trials don't need to be identical
    data1.trialinfo(1,1) = data1.trialinfo(1,1) + 1;
  end
  
  %% remerge each trial file into one
  
  dataout = ft_appenddata([],data1,data2);
  
  %% assert identity to input data
  
  if mismatch
    assert(~isfield(dataout,'trialinfo'),'mismatched trialinfo was not removed')
    disp('mismatched trialinfo was successfully removed upon merge by trials')
  else
    % this time we have to correct for the subtle mismatch that was supposed to
    % pass through undetected
    dataout.trialinfo(1,1) = dataout.trialinfo(1,1) - 1;
    assert(isequaln(data.trialinfo,dataout.trialinfo),'something changed from input to output data')
    % data.hdr not included, data.cfg should change
    disp('trialinfo passed through split & merge by trials unchanged')
  end
end
end