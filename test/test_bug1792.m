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

dataset = '/home/common/matlab/fieldtrip/data/test/bug1792/matteo_cHPI.fif';

hdr = ft_read_header(dataset);

offset   = 3; % seconds

begsample = ([0 1 2 3 4 5]+offset) * hdr.Fs+1;
endsample = begsample + hdr.Fs-1;

cfg = [];
cfg.trl = [begsample(:) endsample(:)];
cfg.trl(:,3) = 0;
cfg.dataset = dataset;
cfg.channel = 'MEGMAG';
data = ft_preprocessing(cfg);

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.foilim = [1 500];
cfg.taper = 'hanning';
freq = ft_freqanalysis(cfg, data);

cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'log10';
freqlog = ft_math(cfg, freq);

cfg = [];
cfg.layout = 'neuromag306mag.lay';
cfg.interactive = 'yes';
ft_multiplotER(cfg, freq);
% ft_multiplotER(cfg, freqlog);


cfg = [];
cfg.dataset = dataset;
cfg.bufferdata = 'first';
cfg.offset = 5;
cfg.jumptoeof = 'no';
cfg.channel = 'MEGMAG';
cfg.coilfreq = [293, 307, 314, 321, 328];
ft_realtime_coillocalizer(cfg);

