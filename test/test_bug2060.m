function test_bug2060

% MEM 1500mb
% WALLTIME 00:03:01

% TEST test_bug2060
% TEST ft_read_event read_neuralynx_nev

filename = dccnfilename('/home/common/matlab/fieldtrip/data/test/bug2060/Events.Nev');

event = ft_read_event(filename);

assert(~isempty(event));

% display
event
event(1)

