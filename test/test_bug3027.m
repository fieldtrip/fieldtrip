function test_bug3027

% WALLTIME 00:10:00
% MEM 2gb

% TEST ft_filetype ft_read_header ft_read_data ft_read_event homer2opto

%%

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug3027/janny/data_janny4juni.txt');
assert(ft_filetype(filename, 'ascii_txt'));
try
  hdr = ft_read_header(filename, 'headerformat', 'bucn_nirs');
  % it should error, the file does not have channel labels
catch
  % this is ok
end

%%

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug3027/janny/being_imitated.sd');
assert(ft_filetype(filename, 'homer_sd'));

%%

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug3027/homer_SampleData/Example1_Simple_Probe/Simple_Probe.nirs');
assert(ft_filetype(filename, 'homer_nirs'));

hdr = ft_read_header(filename);
dat = ft_read_data(filename);

% channel selection
dat1 = ft_read_data(filename, 'chanindx', 1);
dat2 = ft_read_data(filename, 'chanindx', 2);
assert(isequal(dat1, dat(1,:)));
assert(isequal(dat2, dat(2,:)));

% time selection
dat3 = ft_read_data(filename, 'begsample', 33, 'endsample', inf);
assert(isequal(dat3, dat(:,33:end)));
dat4 = ft_read_data(filename, 'begsample', 1, 'endsample', 66);
assert(isequal(dat4, dat(:,1:66)));
dat5 = ft_read_data(filename, 'begsample', 33, 'endsample', 66);
assert(isequal(dat5, dat(:,33:66)));

%%

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug3027/homer_SampleData/Example1_Simple_Probe/Simple_Probe.nirs');

cfg = [];
cfg.dataset = filename;
cfg.viewmode = 'vertical';
cfg.preproc.demean = 'yes';
ft_databrowser(cfg);


%%

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug3027/janny/pp2_13082013.nirs');

hdr = ft_read_header(filename);
dat = ft_read_data(filename);
evt = ft_read_event(filename);

plot([evt.sample], [evt.value], '.')
