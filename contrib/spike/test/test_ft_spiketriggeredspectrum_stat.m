function test_ft_spiketriggeredspectrum_stat()

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_ft_spiketriggeredspectrum_stat
% TEST ft_spiketriggeredspectrum_stat

nSpikes = 10000;
randPhases = [];
kappa = 1;
nChans = 3;
for iChan = 1:nChans
  for iFreq = 1:6
    randPhases(:,iChan,iFreq) = randnwrap(nSpikes+1,kappa/iChan);
  end
end
figure, rose(randPhases(:,1,1))
sts = [];
sts.fourierspctrm{1} = rand(size(randPhases)).*exp(1i*randPhases); % random weighting
sts.trial{1}     = sort(ceil(rand(1,nSpikes)*100));
sts.time{1}         = rand(1,nSpikes);
for iChan = 1:nChans
  sts.lfplabel{iChan} = strcat('chan', num2str(iChan));
end
sts.cfg = [];
sts.dimord         = 'rpt_chan_freq';
sts.freq           = 10:10:60;
sts.label   = {'unit1'};
sts.trialtime      = repmat([0,1],[100 1]);

cfg = [];
cfg.timwin = 0.5;
cfg.winstepsize = 0.001;
cfg.method = 'ppc0';
sts_tfr = ft_spiketriggeredspectrum_stat(cfg,sts);
figure, imagesc(sts_tfr.time, sts_tfr.freq, squeeze(sts_tfr.ppc0(1,:,:))), colorbar
%
cfg.method = 'ppc1';
sts_tfr = ft_spiketriggeredspectrum_stat(cfg,sts);

figure, imagesc(sts_tfr.time, sts_tfr.freq, squeeze(sts_tfr.ppc1(1,:,:))), colorbar
%
cfg = [];
cfg.method = 'ppc1';
sts_tfr = ft_spiketriggeredspectrum_stat(cfg,sts);

figure, plot(sts_tfr.freq, squeeze(sts_tfr.ppc1(1,:,:))), colorbar

%%
% now make a case where no phase locking should be present

% TEST test_ft_spiketriggeredspectrum_tfr
% ft_spiketriggeredspectrum_tfr
nSpikes = 10000;
randPhases = [];
kappa = 1;
nChans = 3;
for iChan = 1:nChans
  for iFreq = 1:6
    randPhases(:,iChan,iFreq) = randnwrap(nSpikes+1,0.001/iChan);
  end
end
figure, rose(randPhases(:,1,1))
sts = [];
sts.fourierspctrm{1} = rand(size(randPhases)).*exp(i*randPhases); % random weighting
sts.trial{1}     = sort(ceil(rand(1,nSpikes)*100));
sts.time{1}         = rand(1,nSpikes);
for iChan = 1:nChans
  sts.lfplabel{iChan} = strcat('chan', num2str(iChan));
end
sts.cfg = [];
sts.dimord         = 'rpt_chan_freq';
sts.freq           = 10:10:60;
sts.label   = {'unit1'};
sts.trialtime      = repmat([0,1],[100 1]);

cfg = [];
cfg.timwin = 0.5;
cfg.winstepsize = 0.001;
cfg.method = 'ppc0';
sts_tfr = ft_spiketriggeredspectrum_stat(cfg,sts);

figure, imagesc(sts_tfr.time, sts_tfr.freq, squeeze(sts_tfr.ppc0(1,:,:))), colorbar
function r=randnwrap(n,k)
r=angle(exp(1i*randn(n,1)*sqrt(1/k)));
