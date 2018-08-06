function test_bug1266

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug_1266
% TEST ft_read_header ft_read_data ft_read_event read_biosig_data read_biosig_header

filename  = fullfile(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1266/sample.gdf'));
filename1 = fullfile(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1266/sample_1.gdf'));

hdr  = ft_read_header(filename);   % this should return the combined header of all three files
hdr1 = ft_read_header(filename1);  % this is only for one file

% the three files are all identical
if hdr.nSamples~=hdr1.nSamples*3
  error('incorrect number of samples')
end

dat  = ft_read_data(filename);
dat1 = ft_read_data(filename1);

if size(dat,2)~=hdr.nSamples
  error('incorrect number of samples')
end

if size(dat1,2)~=hdr1.nSamples
  error('incorrect number of samples')
end

% this is a carefully crafted test that relies on the fileset containing three times the same data
begsample = 97419;  % last sample from 1st file
endsample = 194839; % first sample from 3rd file
dat = ft_read_data(filename, 'begsample', begsample, 'endsample', endsample, 'chanindx', [1 2 3 4 1]);

assert(isequal(dat(1,:), dat(5,:)));
assert(isequal(dat(:,1), dat(:,end-1)));

evt  = ft_read_event(filename);
evt1 = ft_read_event(filename1);
