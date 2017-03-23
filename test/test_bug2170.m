function test_bug2170

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_filetype ft_read_event

% the 3rd dataset is different in the sense that it is a *.fif file accompanied
% by an *.eve file. Those files are also used in the Tristan babysquid74 system.

filename1 = dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/neuromag306/raw.fif');
filename2 = dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/neuromag306/run_01_raw.fif');
filename3 = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2170/an05a4_ss.fif');

assert(ft_filetype(filename1, 'neuromag_fif'));
assert(ft_filetype(filename2, 'neuromag_fif'));
assert(ft_filetype(filename3, 'neuromag_fif'));

event1 = ft_read_event(filename1);
event2 = ft_read_event(filename2);
event3 = ft_read_event(filename3);

assert(length(unique({event1.type}))==8+1); % 'STI001' 'STI002' 'STI003' 'STI004' 'STI009' 'STI010' 'STI011' 'STI101' and an additional 'Trigger'
assert(length(unique({event2.type}))==1);   % 'STI101'
assert(length(unique({event3.type}))==2);   % 'STI101' and 'trigger'

