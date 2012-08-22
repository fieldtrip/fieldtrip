function test_bug1490

% TEST test_bug1490
% TEST ft_read_header ft_read_data loadcnt

datadir       = '/home/common/matlab/fieldtrip/data/test/bug1490';
referencefile = '/home/common/matlab/fieldtrip/data/test/bug1490.mat';

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
filename = fullfile(datadir, filelist{i})

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

%% this is a second test, pertaining to http://bugzilla.fcdonders.nl/show_bug.cgi?id=1490#c11

filename = fullfile(datadir, filelist{2});
hdr = ft_read_header(filename);
% try to titrate the option related to the error, the first ones all fail
dat = ft_read_data(filename, 'header',  hdr, 'begsample',  105704, 'endsample', 206440, 'chanindx', 1:hdr.nChans, 'checkboundary', 1);
dat = ft_read_data(filename, 'header',  hdr, 'begsample',  105704, 'endsample', 206440, 'chanindx', 1:hdr.nChans);
dat = ft_read_data(filename, 'header',  hdr, 'begsample',  105704, 'endsample', 206440, 'chanindx', 1);
dat = ft_read_data(filename,                 'begsample',  105704, 'endsample', 206440, 'chanindx', 1);
dat = ft_read_data(filename,                 'begsample',  105704, 'endsample', inf,    'chanindx', 1); % this failed on 22 aug 2012
dat = ft_read_data(filename,                 'begsample',  105704, 'endsample', inf);                   % this failed on 22 aug 2012
dat = ft_read_data(filename,                 'begsample',       1, 'endsample', inf,    'chanindx', 1); % this worked on 22 aug 2012
dat = ft_read_data(filename,                 'begsample',       1, 'endsample', inf);                   % this worked on 22 aug 2012

% this particular file seems to be in 40-sample blocks
dataformat= 'ns_cnt16'; % ns_cnt16, ns_cnt32
dat = ft_read_data(filename, 'dataformat', dataformat, 'begsample', 1, 'endsample', 80);
assert(size(dat,2)==80, 'incorrect number of samples');
dat = ft_read_data(filename, 'dataformat', dataformat, 'begsample', 1, 'endsample', 40);
assert(size(dat,2)==40, 'incorrect number of samples');
dat = ft_read_data(filename, 'dataformat', dataformat, 'begsample', 1, 'endsample', 39);
assert(size(dat,2)==39, 'incorrect number of samples');
dat = ft_read_data(filename, 'dataformat', dataformat, 'begsample', 1, 'endsample', 41);
assert(size(dat,2)==41, 'incorrect number of samples');
dat = ft_read_data(filename, 'dataformat', dataformat, 'begsample', 2, 'endsample', 41);
assert(size(dat,2)==40, 'incorrect number of samples');



