function test_ft_spiketriggeredspectrum_tfr()

% TEST test_ft_spiketriggeredspectrum_tfr
% ft_spiketriggeredspectrum_tfr
nSpikes = 10000;
randPhases = [];
kappa = 1;
nChans = 3;
for iChan = 1:nChans
  for iFreq = 1:6
    randPhases(:,iChan,iFreq) = generatevonmisesrandvar(nSpikes+1,kappa/iChan);
  end
end
figure, rose(randPhases(:,1,1))
sts = [];
sts.fourierspctrm{1} = rand(size(randPhases)).*exp(i*randPhases); % random weighting
sts.trial{1}     = sort(ceil(rand(1,nSpikes)*100));
sts.time{1}         = rand(1,nSpikes);
for iChan = 1:nChans
  sts.label{iChan} = strcat('chan', num2str(iChan));
end
sts.cfg = [];
sts.dimord         = 'rpt_chan_freq';
sts.freq           = 10:10:60;
sts.spikechannel   = {'unit1'};
sts.trialtime      = [0,1];

cfg = [];
sts_tfr = ft_spiketriggeredspectrum_tfr(cfg,sts);

figure, imagesc(sts_tfr.time, sts_tfr.freq, squeeze(sts_tfr.ppc0(:,1,:))'), colorbar

%% now make a case where no phase locking should be present

clear
nSpikes = 10000;
randPhases = [];
kappa = 1;
nChans = 3;
for iChan = 1:nChans
  for iFreq = 1:6
    randPhases(:,iChan,iFreq) = generatevonmisesrandvar(nSpikes+1,0.0001);
  end
end
sts = [];
sts.fourierspctrm{1} = rand(size(randPhases)).*exp(i*randPhases); % random weighting
sts.trial{1}     = ones(1,nSpikes);
sts.time{1}         = rand(1,nSpikes);
for iChan = 1:nChans
  sts.label{iChan} = strcat('chan', num2str(iChan));
end
sts.cfg = [];
sts.dimord         = 'rpt_chan_freq';
sts.freq           = 10:10:60;
sts.spikechannel   = {'unit1'};
sts.trialtime      = [0,1];
cfg = [];
sts_tfr = ft_spiketriggeredspectrum_tfr(cfg,sts);

figure, imagesc(sts_tfr.time, sts_tfr.freq, squeeze(sts_tfr.ppc0(:,1,:))'), colorbar

