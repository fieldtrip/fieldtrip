function test_pull694

% WALLTIME 00:20:00
% MEM 3gb
% DEPENDENCY read_eyelink_asc

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/original/eyetracker/eyelink'));

filename = {
  'BP2D1_2.asc'
  'ashcal.asc'
  'sub-v31_acq-l2_task-hyperlink_eyetrack.asc'
  'example_file.asc'
  };

% the hard-coded numbers below related to the size of the data were determined prior to the suggested change

%%
i = 1;
hdr = ft_read_header(filename{i});
dat = ft_read_data(filename{i});
evt = ft_read_event(filename{i});
assert(isequal(size(dat), [4 782756]));
assert(numel(evt)==11);

%%
i = 2;
hdr = ft_read_header(filename{i});
dat = ft_read_data(filename{i});
evt = ft_read_event(filename{i});
assert(isequal(size(dat), [7 104517]));
assert(numel(evt)==210);

%%
i = 3;
hdr = ft_read_header(filename{i});
dat = ft_read_data(filename{i});
evt = ft_read_event(filename{i});
assert(isequal(size(dat), [4 1534424]));
assert(numel(evt)==1975);


