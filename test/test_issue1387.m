function test_issue1387

%datadir = dccnpath('/home/common/matlab/fieldtrip/data/');
fileName       = 'Test_4.edf';
bitChannels    = 3:10;
cfg            = [];
cfg.dataset    = fileName;
cfg.header     = ft_read_header(cfg.dataset);
labels         = cfg.header.label(bitChannels);
data_events	   = ft_read_event(cfg.dataset,'header',cfg.header,...
	'detectflank','up','chanindx',bitChannels,'threshold','5*nanmedian');

cfg.continuous = 'yes';
cfg.channel    = 'all';
data           = ft_preprocessing(cfg);

idx            = cellfun(@(x)strcmpi(x,labels{1}),{data_events.type});
event.label	   = labels{1};
event.idx		= find(idx==1);
event.evnt		= data_events(event.idx);
event.samples	= [event.evnt.sample];
event.times	   = data.time{1}(event.samples);

tm             = data.time{1};
ch             = data.trial{1}(3,:);
baseline       = median(ch(1:100));
ch             = ch - baseline;
ch             = ch / max(ch);
figure;
plot(tm,ch,'k-'); 
hold on;
y = repmat(0.5, [1 length(event.times)]);
plot(event.times,y,'r.','MarkerSize',12)
xlim([10 20]);
ylim([-0.05 1.05]);
title('TTL Bit channel 1 should have 78 events')
xlabel('Time (s)')
ylabel('Normalised Amplitude')

assert(length(event.times)==78,'Number of events should be 78!')
