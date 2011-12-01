function test_ft_spike_phaselockstat()

% TEST test_ft_spike_phaselockstat
% ft_spike_phaselockstat

clear
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
  sts.lfplabel{iChan} = strcat('chan', num2str(iChan));
end
sts.cfg = [];
sts.dimord         = 'rpt_chan_freq';
sts.freq           = 10:10:60;
sts.label   = {'unit1'};
sts.trialtime      = [0,1];

cfgTfr = [];
tfr = ft_spike_phaselockstat(cfgTfr,sts);

tfr.ang
tfr.ppc0 % this should give high numbers
tfr.ppc1 % this should give similar values
tfr.ppc2 % this should give similar values
tfr.plv./tfr.ppc1 % should be about sqrt of ppc
tfr.ral % should be significant
if unique(tfr.dofspike)~=10000 % check
  error('dof of spikes incorrect')
end
if unique(tfr.doftrial)~=100 % check
  error('dof of trials incorrect')
end

%%
% now make a case with only a single trial and check that ppc1 becomes nan
% but ppc 0 not

clear
nSpikes = 10000;
randPhases = [];
kappa = 1;
nChans = 3;
for iChan = 1:nChans
  for iFreq = 1:6
    randPhases(:,iChan,iFreq) = generatevonmisesrandvar(nSpikes+1,kappa/iChan);
  end
end
sts = [];
sts.fourierspctrm{1} = rand(size(randPhases)).*exp(i*randPhases); % random weighting
sts.trial{1}     = ones(1,nSpikes);
sts.time{1}         = rand(1,nSpikes);
for iChan = 1:nChans
  sts.lfplabel{iChan} = strcat('chan', num2str(iChan));
end
sts.cfg = [];
sts.dimord         = 'rpt_chan_freq';
sts.freq           = 10:10:60;
sts.label   = {'unit1'};
sts.trialtime      = [0,1];

cfgTfr = [];
tfr = ft_spike_phaselockstat(cfgTfr,sts);

tfr.ang
tfr.ppc0 % this should give high numbers
tfr.ppc1 % this should give similar values
tfr.ppc2 % this should give similar values
tfr.plv./tfr.ppc1 % should be about sqrt of ppc
tfr.ral % should be significant
if unique(tfr.dofspike)~=10000 % check
  error('dof of spikes incorrect')
end
if unique(tfr.doftrial)~=1 % check
  error('dof of trials incorrect')
end

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
  sts.lfplabel{iChan} = strcat('chan', num2str(iChan));
end
sts.cfg = [];
sts.dimord         = 'rpt_chan_freq';
sts.freq           = 10:10:60;
sts.label          = {'unit1'};
sts.trialtime      = [0,1];

cfgTfr = [];
tfr = ft_spike_phaselockstat(cfgTfr,sts);

tfr.ang
tfr.ppc0 % this should give high numbers
tfr.ppc1 % this should give similar values
tfr.ppc2 % this should give similar values
tfr.plv./tfr.ppc1 % should be about sqrt of ppc
tfr.ral % should be significant
if unique(tfr.dofspike)~=10000 % check
  error('dof of spikes incorrect')
end
if unique(tfr.doftrial)~=1 % check
  error('dof of trials incorrect')
end

