function test_bug1792

% MEM 1500mb
% WALLTIME 00:10:00

% this script should not be included in the regression test (yet)
return

% TEST ft_realtime_headlocalizer

fieldtripdir = mfilename('fullpath');
fieldtripdir = fileparts(fieldtripdir); % strip the filename
fieldtripdir = fileparts(fieldtripdir); % strip the test part
disp(fieldtripdir)
addpath(fullfile(fieldtripdir, 'realtime/online_meg'))

dataset = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1792/20130418_test_cHPI.fif');

cfg = [];
cfg.dataset = dataset;
cfg.gradfile = dataset;
cfg.bufferdata = 'first';
cfg.coilfreq = [293, 307, 314, 321];%, 328];
%ft_realtime_coillocalizer(cfg);
ft_realtime_headlocalizer(cfg);

% % stream the testdata to ft buffer in one MATLAB session
% cfg = [];
% cfg.source.datafile = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1792/20130418_test_cHPI.fif');
% cfg.target.datafile = 'buffer://localhost:1972';
% cfg.speed = 1/10;
% ft_realtime_fileproxy(cfg)
% 
% % and then stream from in another matlab session
% cfg = [];
% cfg.dataset = 'buffer://localhost:1972';
% cfg.gradfile = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1792/20130418_test_cHPI.fif');
% cfg.coilfreq = [293, 307, 314, 321];%, 328];
% ft_realtime_headlocalizer(cfg)
