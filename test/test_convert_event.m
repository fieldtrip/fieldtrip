function test_convert_event

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY artifact2boolvec artifact2event artifact2trl boolvec2artifact boolvec2event boolvec2trl event2artifact event2boolvec event2trl trl2artifact trl2boolvec trl2event

% The functionality of converting between representations used to be covered by
% fieldtrip/private/convert_event but has been split over multiple functions to
% provide a better overview of the different conversions.

[ftver, ftpath] = ft_version;
cd(fullfile(ftpath, 'private'));

if false
  % this section helps to get all files opened in the editor
  edit	artifact2boolvec.m
  edit	artifact2event.m
  edit	artifact2trl.m
  edit	boolvec2artifact.m
  edit	boolvec2event.m
  edit	boolvec2trl.m
  edit	event2artifact.m
  edit	event2boolvec.m
  edit	event2trl.m
  edit	trl2artifact.m
  edit	trl2boolvec.m
  edit	trl2event.m
end

% there are 4 representations, hence there are 4x3=12 conversions to be tested

%%

artifact = [
  1 9
  11 19
  21 29
  31 39
  41 49
  ];

if true
  artifact = array2table(artifact);
  artifact.Properties.VariableNames = {'begsample', 'endsample'};
  artifact.type = {'a', 'b', 'c', 'd', 'e'}';
  artifact.value = {1, '2', 3, 4, 5}';
end

boolvec   = artifact2boolvec(artifact, 'endsample', 100)
event     = artifact2event(artifact)
trl       = artifact2trl(artifact)

% convert back
artifact1  = boolvec2artifact(boolvec)
artifact2  = event2artifact(event)
artifact3  = trl2artifact(trl)

%%

boolvec = zeros(1,100);
boolvec(1:10)  = 1;
boolvec(21)    = 2;

artifact  = boolvec2artifact(boolvec)
event     = boolvec2event(boolvec)
trl       = boolvec2trl(boolvec)

% convert back
boolvec1 = artifact2boolvec(artifact, 'endsample', 100)
boolvec2 = event2boolvec(event, 'endsample', 100)
boolvec3 = trl2boolvec(trl, 'endsample', 100)

%%

event = [];

event(1).type     = 'trigger';
event(1).value    = 1;
event(1).sample   = 1;
event(1).duration = [];
event(1).offset   = 0;

event(2).type     = 'trigger';
event(2).value    = 2;
event(2).sample   = 10;
event(2).duration = 2;
event(2).offset   = [];

artifact  = event2artifact(event)
boolvec   = event2boolvec(event, 'endsample', 100)
trl       = event2trl(event)

% convert back
event1 = artifact2event(artifact)
event2 = boolvec2event(boolvec)
event3 = trl2event(trl)

%%

trl = [
  1 9 0
  11 19 0
  21 29 0
  31 39 0
  41 49 0
  ];

if false
  trl = array2table(trl);
  trl.Properties.VariableNames = {'begsample', 'endsample', 'offset'};
  trl.type = {'a', 'b', 'c', 'd', 'e'}';
  trl.value = {1, '2', 3, 4, 5}';
end

artifact  = trl2artifact(trl)
boolvec   = trl2boolvec(trl, 'endsample', 100)
event     = trl2event(trl)

% convert back
trl1 = artifact2trl(artifact)
trl2 = boolvec2trl(boolvec)
trl3 = event2trl(event)
