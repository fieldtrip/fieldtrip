function test_bug1490

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_read_header ft_read_data loadcnt

datadir       = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1490');
referencefile = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1490.mat');

filelist = {
  '0500.cnt'
  'cba1ff01.cnt'
  'dronba4dh.cnt'
  'Subject1_MP.cnt'
  '1pas102_working_memory.cnt'
  'CS14_Sess1_V1_short-block.cnt'
  'sub1E3a.cnt'
  'test.cnt'
  };

% this file contains the reference solution that has been visually checked for correctness
try
  reference = load(referencefile);
  fprintf('SUCCESS: read the reference data from %s\n', referencefile);
catch
  reference = [];
end

hdr = {};
dat = {};

for i=1:length(filelist)
filename = fullfile(datadir, filelist{i});

% read the header and the first 5 seconds of the data
hdr{i} = ft_read_header(filename);
disp(hdr{i});

begsample = 1;
endsample = round(5*hdr{i}.Fs);
dat{i} = ft_read_data(filename, 'begsample', begsample, 'endsample', endsample);

if ~isempty(reference)
% compare to the reference solution
assert(isequal(hdr{i}.Fs,        reference.hdr{i}.Fs));
assert(isequal(hdr{i}.nSamples,  reference.hdr{i}.nSamples));
assert(isequal(hdr{i}.nChans,    reference.hdr{i}.nChans));
assert(isequal(dat{i},           reference.dat{i}));
end
end

if isempty(reference)
  fprintf('WARNING: writing the reference data to %s\n', referencefile);
  save(referencefile, 'hdr', 'dat');
end



