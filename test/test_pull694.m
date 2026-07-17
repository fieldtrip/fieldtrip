function test_pull694

% WALLTIME 00:20:00
% MEM 6gb
% DEPENDENCY read_eyelink_asc
% DATA private

cd(dccnpath('/project/3031000.02/test/original/eyetracker/eyelink'));

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
assert(isequal(size(dat), [5 782756]));
assert(sum(strcmp({evt.type}','input'))==11);

%%
i = 2;
hdr = ft_read_header(filename{i});
dat = ft_read_data(filename{i});
evt = ft_read_event(filename{i});
assert(isequal(size(dat), [11 104517]));
assert(sum(strcmp({evt.type}','input'))==210);

%%
i = 3;
hdr = ft_read_header(filename{i});
dat = ft_read_data(filename{i});
evt = ft_read_event(filename{i});
assert(isequal(size(dat), [5 1534424]));
assert(sum(strcmp({evt.type}','input'))==1975);


