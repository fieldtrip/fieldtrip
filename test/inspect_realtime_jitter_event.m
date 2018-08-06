function inspect_realtime_jitter_event

filename = 'buffer://localhost:1972';

%% write events to the buffer, you should have recording and sine2ft running

% the delay should typically be larger than the blocksize for the data
delay = 0.2;

event = [];
event.type      = 'trigger';
event.sample    = 0;
event.value     = 0;
event.duration  = 0;
event.offset    = 0;

while event.value<50
  event.value = event.value + 1;
  disp(event.value);
  
  ft_write_event(filename, event);
  pause(delay);
end


%% read the recorded events

figure
hold on

% without adjustment
event = ft_read_event('2/0002');
plot(diff([event.sample]./256), 'b.');
std(diff([event.sample]./256))

% with adjustment
event = ft_read_event('3/0001');
plot(diff([event.sample]./256), 'r.');
std(diff([event.sample]./256))
