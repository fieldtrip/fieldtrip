function test_old_ft_freqanalysis

% MEM 1gb
% WALLTIME 00:10:00


function test_ft_freqanalysis
% TEST_FT_FREQANALYSIS
% This script tests the ft_freqanalysis functions using simulated data
% A. Stolk and J.M. Schoffelen

% simulate one second of data, samplefreq = 1200 hz
t = (1:1200)/1200; 
a = cos(2*pi*10*t);
b = sin(2*pi*25*t);
c = a + b;

% simulate preprocessed data
%cfg           = [];
%cfg.layout    = 'CTF275.lay';
%cfg.layout    = prepare_layout(cfg);
%data.label    =  cfg.layout.label(1:273,1);
% data.grad.pnt = zeros(595,3);
% data.grad.ori = zeros(595,3);
% data.grad.tra = zeros(302,595);
data.fsample  = 1200;
data.label    = {'chan01';'cos10';'sin25'};
%for j = 1:273
for j = 1
    %data.trial{1,1}(j,:) = c;
  data.trial{1,1} = [c;a;b];  
  data.time{1,1}(j,:)  = t;
end

% ft_freqnalysis_mtmfft
cfg              = [];
cfg.output       = 'pow';
%cfg.channel      = 'MEG';
cfg.channel      = 'all';
cfg.method       = 'mtmfft';
cfg.taper        = 'hanning';
cfg.foilim       = [1 40];
cfg.keeptrials   = 'no';
cfg.keeptapers   = 'no';
mtmfft           = ft_freqanalysis(cfg, data);

% check whether the powerpeaks are at the given frequencies
if mtmfft.powspctrm(1,9) < mtmfft.powspctrm(1,10) > mtmfft.powspctrm(1,11) && ...
        mtmfft.powspctrm(1,24) < mtmfft.powspctrm(1,25) > mtmfft.powspctrm(1,26); 
else  
    error('test_ft_freqanalysis:notEqual', 'Incorrect output for ft_freqanalysis_mtmfft.');
end

% ft_freqnalysis_mtmconvol
cfg              = [];
cfg.output       = 'pow';
%cfg.channel      = 'MEG';
cfg.channel      = 'all';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.keeptrials   = 'no';
cfg.keeptapers   = 'no';
cfg.foi          = 1:1:40;
cfg.t_ftimwin    = 0.5 * ones(1,length(cfg.foi)); % 500 ms
cfg.toi          = 0.3:0.05:0.75; % center 1/2 second
mtmconvol        = ft_freqanalysis(cfg, data);

% check whether the powerpeaks are at the given frequencies
if mtmconvol.powspctrm(1,9,1) < mtmconvol.powspctrm(1,10,1) > mtmconvol.powspctrm(1,11,1) && ...
        mtmconvol.powspctrm(1,24,1) < mtmconvol.powspctrm(1,25,1) > mtmconvol.powspctrm(1,26,1); 
else  
    error('test_ft_freqanalysis:notEqual', 'Incorrect output for ft_freqanalysis_convol.');
end


% test new implementation specest

% ft_freqnalysis_mtmfft
cfg              = [];
cfg.output       = 'fourier';
cfg.channel      = 'all';
cfg.method       = 'mtmfft';
cfg.taper        = 'hanning';
cfg.foilim       = [1 40];
cfg.keeptrials   = 'yes';
cfg.keeptapers   = 'yes';
mtmfft1          = ft_freqanalysis(cfg, data);
mtmfft2          = ft_freqanalysis(cfg, data, 1);

x = [mtmfft1.fourierspctrm(:,1,10) mtmfft2.fourierspctrm(:,1,10)];
y = [mtmfft1.fourierspctrm(:,1,25) mtmfft2.fourierspctrm(:,1,25)];

% x(2) should only have real component (angle = 0)
% y(2) should only have imag component (negative value, angle = -pi/2)
% abs(x(1))==abs(x(2))
% abs(y(1))==abs(y(2))

% ft_freqnalysis_mtmconvol
cfg              = [];
cfg.output       = 'fourier';
cfg.channel      = 'all';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.keeptrials   = 'no';
cfg.keeptapers   = 'no';
cfg.foi          = 1:1:40;
cfg.t_ftimwin    = 0.5 * ones(1,length(cfg.foi)); % 500 ms
cfg.toi          = 0.25:1./1200:0.75; % center 1/2 second
mtmconvol1        = ft_freqanalysis(cfg, data);
mtmconvol2        = ft_freqanalysis(cfg, data, 1);

x = [squeeze(mtmconvol1.fourierspctrm(:,2,10,:)) ...
     squeeze(mtmconvol2.fourierspctrm(:,2,10,:))];
y = [squeeze(mtmconvol1.fourierspctrm(:,3,25,:)) ...
     squeeze(mtmconvol2.fourierspctrm(:,3,25,:))];

% observations: 
%  Old implementation has 1 Nan at the beginning
%  New implementation has 1 Nan at the end

% expectations:
%  angle(x(1,2)) = -pi, (phase of cosine @10Hz @0.25 s: this is approximately true, but not exact
%  angle(y(1,2)) = 0 (phase of sine @25Hz @0.25 s = 6.25 cycle: this is approximately true, but not exact

% FIXME should we look into this?
% issues are probably related to even numbered t_ftimwins...

% ft_freqnalysis_mtmconvol
cfg              = [];
cfg.output       = 'fourier';
cfg.channel      = 'all';
cfg.method       = 'mtmconvol';
cfg.taper        = 'dpss';
cfg.keeptrials   = 'no';
cfg.keeptapers   = 'no';
cfg.foi          = 1:1:40;
cfg.t_ftimwin    = 0.5 * ones(1,length(cfg.foi)); % 500 ms
cfg.toi          = 0.25:1./1200:0.75; % center 1/2 second
cfg.tapsmofrq    = ones(1,numel(cfg.foi)).*4;
mtmconvol1        = ft_freqanalysis(cfg, data);
mtmconvol2        = ft_freqanalysis(cfg, data, 1);

x = [squeeze(mtmconvol1.fourierspctrm(:,2,10,:)); ...
     squeeze(mtmconvol2.fourierspctrm(:,2,10,:))];
y = [squeeze(mtmconvol1.fourierspctrm(:,3,25,:)); ...
     squeeze(mtmconvol2.fourierspctrm(:,3,25,:))];
