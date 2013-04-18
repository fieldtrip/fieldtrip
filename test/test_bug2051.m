function test_bug2051

% TEST test_bug2051
% TEST ft_math

load(dccnfilename('/home/common/matlab/fieldtrip/data/test/bug2051/source_coh_lft.mat'))

cfg = [];
cfg.parameter = 'avg.pow';
cfg.operation = 'log10';
powlog = ft_math(cfg, source_coh_lft);

% sanity check on some other data
timelock = [];
timelock.pow = randn(1,100).^2;
timelock.label = {'a'};
timelock.time = 1:100;
timelock.dimord = 'chan_time';

<<<<<<< HEAD
cfg = [];
=======
cg = [];
>>>>>>> enhancement - added full provenance to ft_selectdata, implemented support for (new-style) source data in ft_selectdata_new, use rollback_provenance to keep provenance and cfg intact when doing an excursion from ft_math to ft_selectdata, added a 2013x version in ft_datatype_source (to reflect the new-style, still to be discussed with jansch), added a test script for http://bugzilla.fcdonders.nl/show_bug.cgi?id=2053
cfg.parameter = 'pow';
cfg.operation = 'log10';
powlog = ft_math(cfg, timelock);
