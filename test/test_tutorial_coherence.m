function test_tutorial_coherence

% MEM 4500mb
% WALLTIME 00:20:00

% TEST ft_freqanalysis ft_connectivityanalysis ft_multiplotER ft_singleplotER ft_topoplotER ft_sourceanalysis ft_sourceinterpolate ft_prepare_sourcemodel headsurface

addpath(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/coherence'));
addpath(dccnpath('/home/common/matlab/fieldtrip/data/'));

% find the interesting epochs of data
cfg = [];
cfg.trialfun                  = 'trialfun_left';
cfg.dataset                   = dccnpath('/home/common/matlab/fieldtrip/data/SubjectCMC.ds');
cfg = ft_definetrial(cfg);

% detect EOG artifacts in the MEG data
cfg.continuous                = 'yes';
cfg.artfctdef.eog.padding     = 0;
cfg.artfctdef.eog.bpfilter    = 'no';
cfg.artfctdef.eog.detrend     = 'yes';
cfg.artfctdef.eog.hilbert     = 'no';
cfg.artfctdef.eog.rectify     = 'yes';
cfg.artfctdef.eog.cutoff      = 2.5;
cfg.artfctdef.eog.interactive = 'no';
cfg = ft_artifact_eog(cfg);

% detect jump artifacts in the MEG data
cfg.artfctdef.jump.interactive = 'no';
cfg.padding                    = 5;
cfg = ft_artifact_jump(cfg);

% detect muscle artifacts in the MEG data
cfg.artfctdef.muscle.cutoff   = 8;
cfg.artfctdef.muscle.interactive = 'no';
cfg = ft_artifact_muscle(cfg);

% reject the epochs that contain artifacts
cfg.artfctdef.reject          = 'complete';
cfg = ft_rejectartifact(cfg);

% preprocess the MEG data
cfg.demean                    = 'yes';
cfg.dftfilter                 = 'yes';
cfg.channel                   = {'MEG'};
cfg.continuous                = 'yes';
meg = ft_preprocessing(cfg);

% preprocess the EMG data
cfg              = [];
cfg.dataset      = meg.cfg.dataset;
cfg.trl          = meg.cfg.trl;
cfg.continuous   = 'yes';
cfg.demean       = 'yes';
cfg.dftfilter    = 'yes';
cfg.channel      = {'EMGlft' 'EMGrgt'};
cfg.hpfilter     = 'yes';
cfg.hpfreq       = 10;
cfg.rectify      = 'yes';
emg = ft_preprocessing(cfg);

% concatenate the two data-structures into one structure
data = ft_appenddata([], meg, emg);

% visualisation
figure
subplot(2,1,1);
plot(data.time{1},data.trial{1}(77,:));
axis tight;
legend(data.label(77));

subplot(2,1,2);
plot(data.time{1},data.trial{1}(152:153,:));
axis tight;
legend(data.label(152:153));

% spectral analysis: fourier
cfg            = [];
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.foilim     = [5 100];
cfg.tapsmofrq  = 5;
cfg.keeptrials = 'yes';
cfg.channel    = {'MEG' 'EMGlft' 'EMGrgt'};
freqfourier    = ft_freqanalysis(cfg, data);

% spectral analysis: powandcsd
cfg            = [];
cfg.output     = 'powandcsd';
cfg.method     = 'mtmfft';
cfg.foilim     = [5 100];
cfg.tapsmofrq  = 5;
cfg.keeptrials = 'yes';
cfg.channel    = {'MEG' 'EMGlft' 'EMGrgt'};
cfg.channelcmb = {'MEG' 'EMGlft'; 'MEG' 'EMGrgt'};
freq           = ft_freqanalysis(cfg, data);

% compute coherence
cfg            = [];
cfg.method     = 'coh';
cfg.channelcmb = {'MEG' 'EMG'};
fd             = ft_connectivityanalysis(cfg, freq);
fdfourier      = ft_connectivityanalysis(cfg, freqfourier);

% visualisation
cfg                  = [];
cfg.parameter        = 'cohspctrm';
cfg.xlim             = [5 80];
cfg.refchannel       = 'EMGlft';
cfg.layout           = 'CTF151.lay';
cfg.showlabels       = 'yes';
figure; ft_multiplotER(cfg, fd)

cfg.channel = 'MRC21';
figure; ft_singleplotER(cfg, fd);

cfg                  = [];
cfg.parameter        = 'cohspctrm';
cfg.xlim             = [15 20];
cfg.zlim             = [0 0.1];
cfg.refchannel       = 'EMGlft';
cfg.layout           = 'CTF151.lay';
figure; ft_topoplotER(cfg, fd)

% 2 Hz smoothing
cfg            = [];
cfg.output     = 'powandcsd';
cfg.method     = 'mtmfft';
cfg.foilim     = [5 100];
cfg.tapsmofrq  = 2;
cfg.keeptrials = 'yes';
cfg.channel    = {'MEG' 'EMGlft'};
cfg.channelcmb = {'MEG' 'EMGlft'};
freq2          = ft_freqanalysis(cfg,data);

cfg            = [];
cfg.method     = 'coh';
cfg.channelcmb = {'MEG' 'EMG'};
fd2            = ft_connectivityanalysis(cfg,freq2);

cfg               = [];
cfg.parameter     = 'cohspctrm';
cfg.refchannel    = 'EMGlft';
cfg.xlim          = [5 80];
cfg.channel       = 'MRC21';
figure; ft_singleplotER(cfg, fd, fd2);

% 10 Hz smoothing
cfg            = [];
cfg.output     = 'powandcsd';
cfg.method     = 'mtmfft';
cfg.foilim     = [5 100];
cfg.keeptrials = 'yes';
cfg.channel    = {'MEG' 'EMGlft'};
cfg.channelcmb = {'MEG' 'EMGlft'};
cfg.tapsmofrq = 10;
freq10        = ft_freqanalysis(cfg,data);

cfg            = [];
cfg.method     = 'coh';
cfg.channelcmb = {'MEG' 'EMG'};
fd10          = ft_connectivityanalysis(cfg,freq10);

cfg               = [];
cfg.parameter     = 'cohspctrm';
cfg.xlim          = [5 80];
cfg.ylim          = [0 0.2];
cfg.refchannel    = 'EMGlft';
cfg.channel       = 'MRC21';
figure;ft_singleplotER(cfg, fd, fd2, fd10);

% 50 trials
cfg            = [];
cfg.output     = 'powandcsd';
cfg.method     = 'mtmfft';
cfg.foilim     = [5 100];
cfg.tapsmofrq  = 5;
cfg.keeptrials = 'yes';
cfg.channel    = {'MEG' 'EMGlft'};
cfg.channelcmb = {'MEG' 'EMGlft'};
cfg.trials     = 1:50;  
freq50         = ft_freqanalysis(cfg,data);

cfg            = [];
cfg.method     = 'coh';
cfg.channelcmb = {'MEG' 'EMG'};
fd50           = ft_connectivityanalysis(cfg,freq50);

cfg                  = [];
cfg.parameter        = 'cohspctrm';
cfg.xlim             = [5 100];
cfg.ylim             = [0 0.2];
cfg.refchannel       = 'EMGlft';
cfg.channel          = 'MRC21';
figure; ft_singleplotER(cfg, fd, fd50);

% source reconstruction
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'powandcsd';
cfg.foilim     = [18 18];
cfg.tapsmofrq  = 5;
cfg.keeptrials = 'yes';
cfg.channelcmb = {'MEG' 'MEG';'MEG' 'EMGlft'};
freq           = ft_freqanalysis(cfg, data);

cfg           = [];
cfg.method    = 'dics';
cfg.refchan   = 'EMGlft';
cfg.frequency = 18;
cfg.hdmfile   = 'SubjectCMC.hdm';
cfg.inwardshift     = 1;
cfg.grid.resolution = 1;
cfg.grid.unit       = 'cm';
source        = ft_sourceanalysis(cfg, freq);

mri = ft_read_mri('SubjectCMC.mri');
mri = ft_volumereslice([], mri);

cfg            = [];
cfg.parameter  = 'coh';
cfg.downsample = 2;
interp         = ft_sourceinterpolate(cfg, source, mri);

cfg              = [];
cfg.method       = 'ortho';
%cfg.interactive  = 'yes';
cfg.funparameter = 'coh';
figure; ft_sourceplot(cfg, interp);

%--------------------------------
% subfunction
function trl = trialfun_left(cfg)

% read in the triggers and create a trial-matrix
% consisting of 1-second data segments, in which 
% left ECR-muscle is active.

event = ft_read_event(cfg.dataset);
trig  = [event(find(strcmp('backpanel trigger', {event.type}))).value];
indx  = [event(find(strcmp('backpanel trigger', {event.type}))).sample];

%left-condition
sel = [find(trig==1028):find(trig==1029)];

trig = trig(sel);
indx = indx(sel);
 
trl = [];
for j = 1:length(trig)-1
  trg1 = trig(j);
  trg2 = trig(j+1);
  if trg1<=100 && trg2==2080
    trlok      = [[indx(j)+1:1200:indx(j+1)-1200]' [indx(j)+1200:1200:indx(j+1)]'];
    trlok(:,3) = [0:-1200:-1200*(size(trlok,1)-1)]';
    trl        = [trl; trlok];
  end
end
