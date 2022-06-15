function test_issue2048

% WALLTIME 00:10:00
% MEM 2gb

% test data downloaded from here: https://www.teuniz.net/edf_bdf_testfiles/

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/issue2048/test_generator.edf');

%%
% this works

cfg = [];
cfg.channel = 'FP2';
cfg.dataset = filename;
data = ft_preprocessing(cfg);

assert(size(data.trial{1},1)==1)
assert(strcmp(data.label{1}, 'FP2'))

%%
% this also works, i.e., explicitly passing the header

hdr = ft_read_header(filename);
idx_to_load = find(ismember(lower(hdr.label), 'fp2'));
dat         = ft_read_data(filename, 'header', hdr, 'chanindx', idx_to_load);

assert(size(dat,1)==1)
assert(length(hdr.label)==10)

%%
% and this also works, in this case the selection is made early and all channels except for FP2 are hidden from the user

hdr = ft_read_header(filename, 'chanindx', 4);
dat = ft_read_data(filename, 'header', hdr);

assert(size(dat,1)==1)
assert(length(hdr.label)==1)

%%
% but this fails

try
  haserror = false;
  hdr = ft_read_header(filename, 'chanindx', 4);
  dat = ft_read_data(filename, 'header', hdr, 'chanindx', 4);
catch
  haserror = true;
end

% the error is actually correct, since all channels except for FP2 are hidden from the user
% and therefore channel 4 cannot be selected any more; FP2 has become channel 1
assert(haserror)

% this is how it should work

hdr = ft_read_header(filename, 'chanindx', 4);
dat = ft_read_data(filename, 'header', hdr, 'chanindx', 1);

assert(size(dat,1)==1)
assert(length(hdr.label)==1)


%%
% this is where it initially failed according to issue 2048

hdr         = ft_read_header(filename);
idx_to_load = find(ismember(lower(hdr.label), 'fp2')); % this is channel 3 from the selected 10 with 200 Hz (out of 16 channels in the EDF file)
dat         = ft_read_data(filename, 'chanindx', idx_to_load);

assert(length(hdr.label)==10)
assert(size(dat,1)==1)

