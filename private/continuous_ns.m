function [trl] = continuous_ns(cfg)

% CONTINUOUS_NS created a trial definition from a Neuroscan *.cnt file
% which subsequently can be used in the EEG/MEG framework
%
% Use as
%   [trl] = continuous_ns(cfg)
%
% where the configuration should contain 
%   cfg.trialdef.trigger  = number or list with triggers
%   cfg.trialdef.prestim  = pre-stimulus in seconds
%   cfg.trialdef.poststim = post-stimulus in seconds    
%
% See also SINGLETRIAL_NS

% Copyright (C) 2003, Robert Oostenveld

% read the header and the event table from the continuous file
tmp = read_ns_cnt(cfg.datafile, 'ldheaderonly', 1);
fprintf('found %d events in continuous neuroscan file\n', length(tmp.event.frame));
eventindx = find(ismember(tmp.event.stimtype, cfg.trialdef.trigger));
trl(:,1) = tmp.event.frame(eventindx) - round(cfg.trialdef.prestim*tmp.rate) + 1;  % begin sample
trl(:,2) = tmp.event.frame(eventindx) + round(cfg.trialdef.poststim*tmp.rate) + 1; % end sample   
trl(:,3) = -round(cfg.trialdef.prestim*tmp.rate);                                  % offset
fprintf('selected %d events based on triggercode\n', size(trl,1));

return
