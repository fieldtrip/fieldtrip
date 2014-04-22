function test_ft_spike_maketrials()

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_ft_spike_maketrials
% TEST ft_spike_maketrials

spike = [];
tsAll = [];
events = [];
for iTrial = 1:100
  timebeg = 0;
  ts = sort(uint64(timebeg+iTrial*1000000+round(1000000*rand(1,10))));
  tsAll = [tsAll; ts(:)];
  events(iTrial,:) = [timebeg+iTrial*1000000 timebeg+(iTrial+1)*1000000];
end
%%  
spike.timestamp{1} = tsAll;
spike.label{1} = 'spike1';
spike.waveform{1} = rand(32,length(tsAll));
%%
events = uint64(round(events));
trl = events;
%%
cfgC = [];
cfgC.trl = trl;
cfgC.timestampspersecond = 1000000;
cfgC.trl(:,3) = -1*1000000;
cfgC.trl = cfgC.trl(1:50,:);
S = ft_spike_maketrials(cfgC,spike);
