function test_pull1688

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY xdf2fieldtrip fileio/private/sccn_xdf

% description of example data
% contains 5 continuous steams and 156 channels in total
% 128 EEG channels :    1 stream x 128 channels
% 28 motion channels:   4 streams x 7 channels 

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/original/xdf/example_EEG_motion_TUB.xdf');


%% 
% test xdf2fieldtrip : default, import all streams
data     = xdf2fieldtrip(filename);
          
cfg = [];
cfg.viewmode = 'vertical';
cfg.preproc.demean = 'yes';
cfg.ylim = [-20 20];
ft_databrowser(cfg, data); % this should import all 156 continuous channels

%%
% test xdf2fieldtrip : selecting streams by effective sampling rate range
% the EEG stream in the test data has srate of approx. 1000Hz
data     = xdf2fieldtrip(filename,'sraterange',[900 1100]);
          
cfg = [];
cfg.viewmode = 'vertical';
cfg.preproc.demean = 'yes';
cfg.ylim = [-20 20];
ft_databrowser(cfg, data); % this should display the 128 eeg channels

%%
% test xdf2fieldtrip : selecting multiple streams by a single shared keyword in stream names
% this will import both "LeftFoot" and "RightFoot" streams containing the keyword "Foot"
data    = xdf2fieldtrip(filename,'streamkeywords',{'Foot'}); 

cfg = [];
cfg.viewmode = 'vertical';
cfg.preproc.demean = 'yes';
cfg.ylim = [-20 20];
ft_databrowser(cfg, data); % this should display 14 motion channels from 2 feet streams


%%
% test xdf2fieldtrip : selecting multiple streams by multiple keywords in stream names 
% this will import all 4 motion streams
data = xdf2fieldtrip(filename,'streamkeywords',{'PlayerTransform', 'LeftFoot', 'RightFoot', 'Torso'});

cfg = [];
cfg.viewmode = 'vertical';
cfg.preproc.demean = 'yes';
cfg.ylim = [-20 20];
ft_databrowser(cfg, data); % this should display 28 motion channels from all 4 motion streams

end
