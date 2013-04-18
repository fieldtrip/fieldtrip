function test_bug2086

% TEST test_bug2086
% TEST ft_databrowser

warning('this test should not run automatically');
return

load(dccnfilename('/home/common/matlab/fieldtrip/data/test/bug2086.mat'));

cfg = [];
<<<<<<< HEAD
=======
%cfg.dataset = sprintf('%s', 'data/' , subj{isub} ,'/', cond{icond},  set{1},time{itime}  ,'_all.mat');
>>>>>>> enhancement - added test script for issue #2086
cfg.continuous = 'no';
cfg.preproc.channel = {'all'  };
cfg.preproc.demean = 'yes';
cfg.viewmode = 'butterfly';
cfg.ylim = [-200 200];
<<<<<<< HEAD
=======
%cfg.layout = 'J:\Projects\ANT_trial\Export\customANT.lay';     
>>>>>>> enhancement - added test script for issue #2086
cfg.channel = 'all';
cfg = ft_databrowser(cfg,data);
