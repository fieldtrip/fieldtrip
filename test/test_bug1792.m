function test_bug1792

% this script should not be included in the regression test (yet)
return

% TEST test_bug1792
% TEST ft_realtime_headlocalizer

fieldtripdir = mfilename('fullpath');
fieldtripdir = fileparts(fieldtripdir); % strip the filename
fieldtripdir = fileparts(fieldtripdir); % strip the test part
disp(fieldtripdir)
addpath(fullfile(fieldtripdir, 'realtime/online_meg'))

dataset = '/home/common/matlab/fieldtrip/data/test/bug1792/20130418_test_cHPI.fif';

cfg = [];
cfg.dataset = dataset;
cfg.bufferdata = 'first';
cfg.coilfreq = [293, 307, 314, 321];%, 328];
%ft_realtime_coillocalizer(cfg);
ft_realtime_headlocalizer(cfg);

