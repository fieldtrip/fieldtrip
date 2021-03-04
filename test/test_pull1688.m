function test_pull1688

% WALLTIME 00:01:00
% MEM 1gb
% DEPENDENCY xdf2fieldtrip fileio/private/sccn_xdf

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/original/xdf/example_EEG_motion_TUB.xdf');


%%
% test xdf2fieldtrip : selecting streams by effective sampling rate range
% the test file contains an EEG stream of srate approx. 1000Hz

eeg     = xdf2fieldtrip(filename,'sraterange',[900 1100]);
          

%%
% test xdf2fieldtrip : selecting multiple streams by a single shared keyword in stream names

feet    = xdf2fieldtrip(filename,'streamkeywords',{'Foot'}); % this will select both LeftFoot and RightFoot streams

%%
% test xdf2fieldtrip : selecting multiple streams by multiple keywords in stream names 

motion = xdf2fieldtrip(filename,'streamkeywords',{'PlayerTransform', 'LeftFoot', 'RightFoot', 'Torso'});


end
