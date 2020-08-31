function test_issue1387

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_read_event ft_read_header ft_preprocessing

%%
% load data and parse the events using '5*nanmedian' threshold on 8 channels
fileName       = dccnpath('/home/common/matlab/fieldtrip/data/test/issue1387/Test_4.edf');
bitChannels    = 3:10;
header         = ft_read_header(fileName);
labels         = header.label(bitChannels);
data_events	   = ft_read_event(fileName,'header',header,'detectflank','up','chanindx',bitChannels,'threshold','5*nanmedian');

%%
% load the raw data
cfg            = [];
cfg.dataset    = fileName;
cfg.continuous = 'yes';
cfg.channel    = 'all';
data           = ft_preprocessing(cfg);

%%
% parse out the event information for the 1st bit channel
idx            = cellfun(@(x)strcmpi(x,labels{1}),{data_events.type});
event.label	   = labels{1};
event.idx      = find(idx==1);
event.evnt     = data_events(event.idx);
event.samples  = [event.evnt.sample];
event.times	   = data.time{1}(event.samples);

%% 
% ASSERT: There were 78 TTLs on the first channel, make sure this is the case, error otherwise.
assert(length(event.times)==78,'Number of events should be 78!')

%%
% Plot the events on top of the raw data for visual inspection
tm             = data.time{1};
ch             = data.trial{1}(3,:);
baseline       = median(ch(1:100));
ch             = ch - baseline;
ch             = ch / max(ch);
figure;
plot(tm,ch,'k-'); 
hold on;
y = repmat(0.5, [1 length(event.times)]);
plot(event.times,y,'ro','MarkerSize',6)
xlim([10 50]);
ylim([-0.05 1.05]);
title('TTL channel 1 should have 78 event markers aligned to the upstates')
xlabel('Time (s)')
ylabel('Normalised Amplitude')

