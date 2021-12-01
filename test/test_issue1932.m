function test_issue1932

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY fieldtrip2homer event2boolvec

load(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1932/data.mat'));
tmp = load(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1932/Events.mat'));
event = tmp.events;

%%

clear values
for i=1:length(event)
  begsample = event(i).sample + event(i).offset;
  endsample = event(i).sample + event(i).offset + event(i).duration;
  values(begsample:endsample) = find(strcmp(event(i).value, unique({event.value})));
end

figure
hold on
plot(1*(values==1), '.')
plot(2*(values==2), '.')
plot(3*(values==3), '.')

axis([56000 76000 0.5 3.5]);
legend(unique({event.value}))

%%

nirs = fieldtrip2homer(Selection, 'event', event);

figure
imagesc(nirs.s)

% there should not be any overlap in events
assert(~any(nirs.s(:,1) & nirs.s(:,2)));
assert(~any(nirs.s(:,1) & nirs.s(:,3)));
assert(~any(nirs.s(:,2) & nirs.s(:,3)));
