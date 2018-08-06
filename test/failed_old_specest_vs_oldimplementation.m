function failed_old_specest_vs_oldimplementation

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_old_specest_vs_oldimplementation

% read in data from /home/common/testdata (I just my own testdata)
% data = .....

%% COHERENCE specest_mtmfft
cfg = [];
cfg.keeptrials = 'yes';
cfg.foilim     = [2 100];
cfg.output     = 'powandcsd';
cfg.taper      = 'dpss';
cfg.tapsmofrq  = 4;
cfg.pad        = 'maxperlen';
cfg.calcdof    = 'yes';
cfg.keeptrials = 'no';
% old style
cfg.method     = 'mtmfft_old';
tic; freqold = ft_freqanalysis(cfg,data); toc
% new style
cfg.method     = 'mtmfft';
tic; freqnew = ft_freqanalysis(cfg,data); toc

% calculate coherence
cfg = [];
cfg.method = 'coh';
connold = ft_connectivityanalysis(cfg,freqold);
connnew = ft_connectivityanalysis(cfg,freqnew);

% plot 'a' chancomb
figure
subplot(3,1,1)
plot(connold.freq,squeeze(connold.cohspctrm(10,:)))
axis tight
xlabel('frequency (Hz)')
ylabel('coh')
title('old implementation')
subplot(3,1,2)
plot(connnew.freq,squeeze(connnew.cohspctrm(10,:)))
axis tight
xlabel('frequency (Hz)')
ylabel('coh')
title('new implementation')
subplot(3,1,3)
hold on
plot(connold.freq,squeeze(connold.cohspctrm(10,:)))
plot(connnew.freq,squeeze(connnew.cohspctrm(10,:)),'r')
legend
axis tight
xlabel('frequency (Hz)')
ylabel('coh')
title('both overlayed, using freq axis of new implementation')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% POWER specest_mtmfft
cfg = [];
cfg.trials     = 1:3;
cfg.keeptrials = 'yes';
cfg.foilim     = [0 256/2];
cfg.output     = 'pow';
cfg.taper      = 'dpss';
cfg.tapsmofrq  = 4;
cfg.pad        = 'maxperlen';
cfg.calcdof    = 'yes';
% old style
cfg.method     = 'mtmfft_old';
tic; freqold = ft_freqanalysis(cfg,data); toc
% new style
cfg.method     = 'mtmfft';
tic; freqnew = ft_freqanalysis(cfg,data); toc

% plot first trial first channel, single timepoint
figure
subplot(3,1,1)
plot(freqold.freq,squeeze(freqold.powspctrm(1,1,:)))
axis tight
xlabel('frequency (Hz)')
ylabel('power')
title('old implementation')
subplot(3,1,2)
plot(freqnew.freq,squeeze(freqnew.powspctrm(1,1,:)))
axis tight
xlabel('frequency (Hz)')
ylabel('power')
title('new implementation')
subplot(3,1,3)
hold on
plot(freqnew.freq,squeeze(freqold.powspctrm(1,1,:)))
plot(freqnew.freq,squeeze(freqnew.powspctrm(1,1,:)),'r')
legend
axis tight
xlabel('frequency (Hz)')
ylabel('power')
title('both overlayed, using freq axis of new implementation')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% COHERENCE specest_mtmconvol
cfg = [];
cfg.trials     = 1:2;
cfg.keeptrials = 'no';
cfg.keeptapers = 'no';
cfg.channel    = 'all';
cfg.channelcmb = {'all' 'all'};
%cfg.pad        = 4.24;
cfg.foi        = 2:2:20;
cfg.toi        = data.time{1}(1):0.050:data.time{1}(end);
cfg.output     = 'powandcsd';
cfg.taper      = 'hanning';
cfg.t_ftimwin  = 4 ./ cfg.foi; %ones(length(cfg.foi),1) .* 1.5;%
cfg.tapsmofrq  = ones(length(cfg.foi),1) * 8;
cfg.calcdof    = 'yes';
cfg.correctt_ftimwin = 'no';
% old style
cfg.method     = 'mtmconvol_old';
tic; freqold = ft_freqanalysis(cfg,data); toc
% new style
cfg.method     = 'mtmconvol';
tic; freqnew = ft_freqanalysis(cfg,data); toc

% calculate coherence
cfg = [];
cfg.method = 'coh';
connold = ft_connectivityanalysis(cfg,freqold);
connnew = ft_connectivityanalysis(cfg,freqnew);



% plot 'a' chancomb at a single timepoint
timepoint = 20;
figure
subplot(3,1,1)
plot(connold.freq,squeeze(connold.cohspctrm(10,:,timepoint)))
axis tight
xlabel('frequency (Hz)')
ylabel('coh')
title('old implementation')
subplot(3,1,2)
plot(connnew.freq,squeeze(connnew.cohspctrm(10,:,timepoint)))
axis tight
xlabel('frequency (Hz)')
ylabel('coh')
title('new implementation')
subplot(3,1,3)
hold on
plot(connold.freq,squeeze(connold.cohspctrm(10,:,timepoint)))
plot(connnew.freq,squeeze(connnew.cohspctrm(10,:,timepoint)),'r')
legend
axis tight
xlabel('frequency (Hz)')
ylabel('coh')
title('both overlayed, using freq axis of new implementation')



% plot 'a' chancomb at a all timepoints
figure
subplot(3,1,1)
imagesc(connold.time,connold.freq,squeeze(connold.cohspctrm(10,:,:)))
axis tight
ylabel('frequency (Hz)')
xlabel('time (s)')
title('old implementation')
colorbar; axis xy; caxis([0 1])
subplot(3,1,2)
imagesc(connnew.time,connnew.freq,squeeze(connnew.cohspctrm(10,:,:)))
axis tight
ylabel('frequency (Hz)')
xlabel('time (s)')
title('new implementation')
colorbar; axis xy; caxis([0 1]);
subplot(3,1,3)
hold on
imagesc(connold.time,connold.freq,squeeze(connold.cohspctrm(10,:,:))-squeeze(connnew.cohspctrm(10,:,:)))
legend
axis tight
ylabel('frequency (Hz)')
xlabel('time (s)')
title('old minus new')
colorbar; caxis([0 1]); axis xy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%% POWER specest_mtmconvol
cfg = [];
cfg.implementation = 'old';
cfg.trials     = 1:3;
cfg.keeptrials = 'yes';
cfg.keeptapers = 'no';
%cfg.pad        = 4.24;
cfg.foi        = 0:2:128;
cfg.toi        = data.time{1}(1):0.0050:data.time{1}(end);
cfg.output     = 'pow';
cfg.taper      = 'hanning';
cfg.t_ftimwin  = 4 ./ cfg.foi;%ones(length(cfg.foi),1) .* 1.5;%
cfg.tapsmofrq  = ones(length(cfg.foi),1) * 8;
cfg.calcdof    = 'yes';
cfg.correctt_ftimwin = 'no';
% new style
cfg.method     = 'mtmconvol';
tic; freqnew = ft_freqanalysis(cfg,data); toc
% old style
cfg.method     = 'mtmconvol_old';
cfg.foi = freqnew.freq;
cfg.t_ftimwin  = 4 ./ cfg.foi;%ones(length(cfg.foi),1) .* 1.5;%
tic; freqold = ft_freqanalysis(cfg,data); toc





% plot first trial first channel, single timepoint
timepoint = 20;
figure
subplot(3,1,1)
plot(freqold.freq,squeeze(freqold.powspctrm(1,1,:,timepoint)))
axis tight
xlabel('frequency (Hz)')
ylabel('power')
title('old implementation')
subplot(3,1,2)
plot(freqnew.freq,squeeze(freqnew.powspctrm(1,1,:,timepoint)))
axis tight
xlabel('frequency (Hz)')
ylabel('power')
title('new implementation')
subplot(3,1,3)
hold on
plot(freqnew.freq,squeeze(freqold.powspctrm(1,1,:,timepoint)))
plot(freqnew.freq,squeeze(freqnew.powspctrm(1,1,:,timepoint)),'r')
legend
axis tight
xlabel('frequency (Hz)')
ylabel('power')
title('both overlayed, using freq axis of new implementation')



% plot first trial first channel, all timepoints
figure
subplot(3,1,1)
imagesc(freqold.time,freqold.freq,squeeze(freqold.powspctrm(1,1,:,:)))
axis tight; axis xy; colorbar
ylabel('frequency (Hz)')
xlabel('time (s)')
title('old implementation')
subplot(3,1,2)
imagesc(freqold.time,freqnew.freq,squeeze(freqnew.powspctrm(1,1,:,:)))
axis tight; axis xy; colorbar
ylabel('frequency (Hz)')
xlabel('time (s)')
title('new implementation')
subplot(3,1,3)
hold on
imagesc(freqold.time,freqold.freq,squeeze(freqold.powspctrm(1,1,:,:)) - squeeze(freqnew.powspctrm(1,1,:,:)))
legend
axis tight; axis xy; colorbar
ylabel('frequency (Hz)')
xlabel('time (s)')
title('old minus new')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%% POWER specest_wavelet
cfg = [];
cfg.trials     = 1;
cfg.keeptrials = 'yes';
cfg.keeptapers = 'no';
%cfg.pad        = 4.24;
cfg.foi        = 2:1:40;
cfg.toi        = 0:0.05:4;%data.time{1}(1):(1/data.fsample):data.time{1}(end);
cfg.output     = 'pow';
cfg.calcdof    = 'yes';
cfg.correctt_ftimwin = 'no';
% new style
cfg.method     = 'wavelet';
tic; freqnew = ft_freqanalysis(cfg,data); toc
% old style
cfg.foi = freqnew.freq;
cfg.method     = 'wltconvol_old';
tic; freqold = ft_freqanalysis(cfg,data); toc



% plot first trial first channel, single timepoint
timepoint = 40;
figure
subplot(3,1,1)
plot(freqold.freq,squeeze(freqold.powspctrm(1,1,:,timepoint)))
axis tight
xlabel('frequency (Hz)')
ylabel('power')
title('old implementation')
subplot(3,1,2)
plot(freqnew.freq,squeeze(freqnew.powspctrm(1,1,:,timepoint)))
axis tight
xlabel('frequency (Hz)')
ylabel('power')
title('new implementation')
subplot(3,1,3)
hold on
plot(freqnew.freq,squeeze(freqold.powspctrm(1,1,:,timepoint)))
plot(freqnew.freq,squeeze(freqnew.powspctrm(1,1,:,timepoint)),'r')
legend
axis tight
xlabel('frequency (Hz)')
ylabel('power')
title('both overlayed, using freq axis of new implementation')



% plot first trial first channel, all timepoints
figure
subplot(3,1,1)
imagesc(freqold.time,freqold.freq,squeeze(freqold.powspctrm(1,1,:,:)))
axis tight; axis xy; colorbar
ylabel('frequency (Hz)')
xlabel('time (s)')
title('old implementation')
subplot(3,1,2)
imagesc(freqold.time,freqnew.freq,squeeze(freqnew.powspctrm(1,1,:,:)))
axis tight; axis xy; colorbar
ylabel('frequency (Hz)')
xlabel('time (s)')
title('new implementation')
subplot(3,1,3)
hold on
imagesc(freqold.time,freqold.freq,squeeze(freqold.powspctrm(1,1,:,:)) - squeeze(freqnew.powspctrm(1,1,:,:)))
legend
axis tight; axis xy; colorbar
ylabel('frequency (Hz)')
xlabel('time (s)')
title('old minus new')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% CSD specest_wavelet
cfg = [];
cfg.trials     = 'all';
cfg.keeptrials = 'yes';
cfg.keeptapers = 'no';
%cfg.pad        = 4.24;
cfg.foi        = 2:2:40;
cfg.toi        = data.time{1}(1):(1/data.fsample):data.time{1}(end);
cfg.output     = 'powandcsd';
cfg.calcdof    = 'yes';
cfg.correctt_ftimwin = 'no';
% new style
cfg.method     = 'wavelet';
tic; freqnew = ft_freqanalysis(cfg,data); toc
% old style
cfg.foi = freqnew.freq;
cfg.method     = 'wltconvol_old';
tic; freqold = ft_freqanalysis(cfg,data); toc




% plot first trial first channel, single timepoint
timepoint = 200;
figure
subplot(3,1,1)
plot(freqold.freq,abs(squeeze(freqold.crsspctrm(1,1,:,timepoint))))
axis tight
xlabel('frequency (Hz)')
ylabel('')
title('old implementation')
subplot(3,1,2)
plot(freqnew.freq,abs(squeeze(freqnew.crsspctrm(1,1,:,timepoint))))
axis tight
xlabel('frequency (Hz)')
ylabel('')
title('new implementation')
subplot(3,1,3)
hold on
plot(freqnew.freq,abs(squeeze(freqold.crsspctrm(1,1,:,timepoint))))
plot(freqnew.freq,abs(squeeze(freqnew.crsspctrm(1,1,:,timepoint))),'r')
legend
axis tight
xlabel('frequency (Hz)')
ylabel('')
title('both overlayed, using freq axis of new implementation')



% plot first trial first channel, all timepoints
figure
subplot(3,1,1)
imagesc(freqold.time,freqold.freq,abs(squeeze(freqold.crsspctrm(1,1,:,:))))
axis tight; axis xy; colorbar
ylabel('frequency (Hz)')
xlabel('time (s)')
title('old implementation')
subplot(3,1,2)
imagesc(freqold.time,freqnew.freq,abs(squeeze(freqnew.crsspctrm(1,1,:,:))))
axis tight; axis xy; colorbar
ylabel('frequency (Hz)')
xlabel('time (s)')
title('new implementation')
subplot(3,1,3)
hold on
imagesc(freqold.time,freqold.freq,abs(squeeze(freqold.crsspctrm(1,1,:,:)) - squeeze(freqnew.crsspctrm(1,1,:,:))))
legend
axis tight; axis xy; colorbar
ylabel('frequency (Hz)')
xlabel('time (s)')
title('old minus new')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% COHERENCE specest_wavelet
cfg = [];
cfg.trials     = 'all';
cfg.keeptrials = 'yes';
cfg.keeptapers = 'no';
cfg.channel    = 'all';
cfg.channelcmb = {'all' 'all'};
%cfg.pad        = 4.24;
cfg.foi        = 2:2:20;
cfg.toi        = data.time{1}(1):0.050:data.time{1}(end);
cfg.output     = 'powandcsd';
cfg.calcdof    = 'yes';
% new style
cfg.method     = 'wavelet';
tic; freqnew = ft_freqanalysis(cfg,data); toc
% old style
cfg.foi = freqnew.freq;
cfg.method     = 'wltconvol_old';
tic; freqold = ft_freqanalysis(cfg,data); toc

% calculate coherence
cfg = [];
cfg.method = 'coh';
connold = ft_connectivityanalysis(cfg,freqold);
connnew = ft_connectivityanalysis(cfg,freqnew);



% plot 'a' chancomb at a single timepoint
timepoint = 20;
figure
subplot(3,1,1)
plot(connold.freq,squeeze(connold.cohspctrm(10,:,timepoint)))
axis tight
xlabel('frequency (Hz)')
ylabel('coh')
title('old implementation')
subplot(3,1,2)
plot(connnew.freq,squeeze(connnew.cohspctrm(10,:,timepoint)))
axis tight
xlabel('frequency (Hz)')
ylabel('coh')
title('new implementation')
subplot(3,1,3)
hold on
plot(connold.freq,squeeze(connold.cohspctrm(10,:,timepoint)))
plot(connnew.freq,squeeze(connnew.cohspctrm(10,:,timepoint)),'r')
legend
axis tight
xlabel('frequency (Hz)')
ylabel('coh')
title('both overlayed, using freq axis of new implementation')



% plot 'a' chancomb at a all timepoints
figure
subplot(3,1,1)
imagesc(connold.time,connold.freq,squeeze(connold.cohspctrm(10,:,:)))
axis tight
ylabel('frequency (Hz)')
xlabel('time (s)')
title('old implementation')
colorbar; axis xy; caxis([0 1])
subplot(3,1,2)
imagesc(connnew.time,connnew.freq,squeeze(connnew.cohspctrm(10,:,:)))
axis tight
ylabel('frequency (Hz)')
xlabel('time (s)')
title('new implementation')
colorbar; axis xy; caxis([0 1]);
subplot(3,1,3)
hold on
imagesc(connold.time,connold.freq,squeeze(connold.cohspctrm(10,:,:))-squeeze(connnew.cohspctrm(10,:,:)))
legend
axis tight
ylabel('frequency (Hz)')
xlabel('time (s)')
title('old minus new')
colorbar; caxis([0 1]); axis xy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% POWER specest_hilbert
cfg = [];
cfg.trials     = 1:3;
cfg.keeptrials = 'yes';
cfg.keeptapers = 'no';
cfg.width      = 1;
cfg.pad        = 'maxperlen';
cfg.foi        = 10:2:40;
cfg.toi        = 0:0.05:4;%data.time{1}(1):(1/data.fsample):data.time{1}(end);
cfg.output     = 'pow';
cfg.calcdof    = 'yes';
cfg.correctt_ftimwin = 'no';
% new style
cfg.method     = 'hilbert';
cfg.filtorder  = round((3 ./ cfg.foi) .* data.fsample);
cfg.filttype   = 'fir';
cfg.width      = 5;

tic; freqnew = ft_freqanalysis(cfg,data); toc
% old style
cfg.foi = freqnew.freq;
cfg.method     = 'mtmconvol';
cfg.taper      = 'hanning';
cfg.t_ftimwin  = 3 ./ cfg.foi;%ones(length(cfg.foi),1) .* 1.5;%
tic; freqold = ft_freqanalysis(cfg,data); toc



% plot first trial first channel, single timepoint
timepoint = 40;
figure
subplot(3,1,1)
plot(freqold.freq,squeeze(freqold.powspctrm(1,1,:,timepoint)))
axis tight
xlabel('frequency (Hz)')
ylabel('power')
title('mtmconvol')
subplot(3,1,2)
plot(freqnew.freq,squeeze(freqnew.powspctrm(1,1,:,timepoint)))
axis tight
xlabel('frequency (Hz)')
ylabel('power')
title('hilbert')
subplot(3,1,3)
hold on
plot(freqnew.freq,squeeze(freqold.powspctrm(1,1,:,timepoint)))
plot(freqnew.freq,squeeze(freqnew.powspctrm(1,1,:,timepoint)),'r')
legend
axis tight
xlabel('frequency (Hz)')
ylabel('power')
title('both overlayed, using freq axis of hilbert')



% plot first trial first channel, all timepoints
figure
subplot(3,1,1)
imagesc(freqold.time,freqold.freq,squeeze(freqold.powspctrm(1,1,:,:)))
axis tight; axis xy; colorbar
ylabel('frequency (Hz)')
xlabel('time (s)')
title('mtmconvol')
subplot(3,1,2)
imagesc(freqold.time,freqnew.freq,squeeze(freqnew.powspctrm(1,1,:,:)))
axis tight; axis xy; colorbar
ylabel('frequency (Hz)')
xlabel('time (s)')
title('hilbert')
subplot(3,1,3)
hold on
imagesc(freqold.time,freqold.freq,squeeze(freqold.powspctrm(1,1,:,:)) - squeeze(freqnew.powspctrm(1,1,:,:)))
legend
axis tight; axis xy; colorbar
ylabel('frequency (Hz)')
xlabel('time (s)')
title('old minus new')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% FOURIER specest_hilbert
cfg = [];
cfg.trials     = 1:3;
cfg.keeptrials = 'yes';
cfg.keeptapers = 'yes';
cfg.width      = 1;
cfg.pad        ='maxperlen';
cfg.foi        = 10:2:40;
cfg.toi        = 0:0.05:4;%data.time{1}(1):(1/data.fsample):data.time{1}(end);
cfg.output     = 'fourier';
cfg.calcdof    = 'yes';
cfg.correctt_ftimwin = 'no';
% new style
cfg.method     = 'hilbert';
cfg.filtorder  = round((3 ./ cfg.foi) .* data.fsample);
cfg.filttype   = 'fir';
cfg.width      = 1;
tic; freqnew = ft_freqanalysis(cfg,data); toc
% old style
cfg.foi = freqnew.freq;
cfg.method     = 'mtmconvol';
cfg.taper      = 'hanning';
cfg.t_ftimwin  = 3 ./ cfg.foi;%ones(length(cfg.foi),1) .* 1.5;%
tic; freqold = ft_freqanalysis(cfg,data); toc



% plot first trial first channel, real part
freq = 3;
figure
subplot(3,1,1)
plot(freqold.time,real(squeeze(freqold.fourierspctrm(1,1,freq,:))))
axis tight
xlabel('frequency (Hz)')
ylabel('amplitude')
title('mtmconvol')
subplot(3,1,2)
plot(freqnew.time,real(squeeze(freqnew.fourierspctrm(1,1,freq,:))))
axis tight
xlabel('frequency (Hz)')
ylabel('amplitude')
title('hilbert')
subplot(3,1,3)
hold on
plot(freqnew.time,real(squeeze(freqold.fourierspctrm(1,1,freq,:))))
plot(freqnew.time,real(squeeze(freqnew.fourierspctrm(1,1,freq,:))),'r')
legend
axis tight
xlabel('frequency (Hz)')
ylabel('amplitude')
title('both overlayed, using freq axis of hilbert')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% FOURIER specest_wavelet VS specest_tfr
cfg = [];
cfg.trials     = 1;
cfg.keeptrials = 'yes';
cfg.keeptapers = 'yes';
cfg.pad        = 6;
cfg.foi        = 2:5:20;
cfg.toi        = 0:0.05:4;%data.time{1}(1):(1/data.fsample):data.time{1}(end);
cfg.output     = 'fourier';
cfg.calcdof    = 'yes';
cfg.correctt_ftimwin = 'no';
% new style
cfg.method     = 'wavelet';
tic; freqwav = ft_freqanalysis(cfg,data); toc
% old style
cfg.foi = freqwav.freq;
cfg.method     = 'tfr';
tic; freqtfr = ft_freqanalysis(cfg,data); toc



% plot first trial first channel, single timepoint
timepoint = 50;
figure
subplot(4,1,1)
plot(freqtfr.freq,real(squeeze(freqtfr.fourierspctrm(1,1,:,timepoint))))
axis tight
xlabel('frequency (Hz)')
ylabel('amplitude')
title('tfr implementation')
subplot(4,1,2)
plot(freqwav.freq,real(squeeze(freqwav.fourierspctrm(1,1,:,timepoint))))
axis tight
xlabel('frequency (Hz)')
ylabel('amplitude')
title('wavelet implementation')
subplot(4,1,3)
plot(freqtfr.freq,angle(squeeze(freqtfr.fourierspctrm(1,1,:,timepoint))))
axis tight
xlabel('frequency (Hz)')
ylabel('phase')
title('tfr implementation')
subplot(4,1,4)
plot(freqwav.freq,angle(squeeze(freqwav.fourierspctrm(1,1,:,timepoint))))
axis tight
xlabel('frequency (Hz)')
ylabel('phase')
title('wavelet implementation')








