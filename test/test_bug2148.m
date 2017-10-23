function test_bug2148

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_connectivitysimulation ft_freqanalysis ft_connectivityanalysis
% TEST ft_connectivityplot ft_freqdescriptives ft_checkdata

% ft_checkdata (ie fixdimord) always wants to reverts dimord to 'chan'
% although 'labelcmb' is present. The only way to make e.g.
% ft_connectivitplot work atm is to use:
% ft_checkdata(data, 'cmbrepresentation', 'full');
% This will make dimord chan_chan_freq and convert labelcmb to label.

cfg             = [];
cfg.ntrials     = 500;
cfg.triallength = 1;
cfg.fsample     = 200;
cfg.nsignal     = 4;
cfg.method      = 'ar';
cfg.params(:,:,1) = [.8   0   0 .2;
  0  .9  .5  0;
  .4   0  .5 .3;
  0  .2   0 .7];
cfg.params(:,:,2) = [-.5   0   0  0;
  0 -.8   0  0;
  0   0 -.2  0;
  0   0   0 .1];
cfg.noisecov      = [.3 0  0  0;
  0 1  0  0;
  0 0 .2  0;
  0 0  0  .4];

data            = ft_connectivitysimulation(cfg);
data.label      = {'AFz', 'Cz', 'Fz', 'Oz'};

% freqanalysis
cfgf           = [];
cfgf.method    = 'mtmfft';
cfgf.output    = 'fourier';
cfgf.tapsmofrq = 2;
cfgf.keeptrials= 'yes';
freq           = ft_freqanalysis(cfgf, data);
ft_freqdescriptives([], freq);

% connectivityanalysis
cfgc           = [];
cfgc.channelcmb = {'all' freq.label{1}};
cfgc.method    = 'coh';
c1             = ft_connectivityanalysis(cfgc, freq);
ft_connectivityplot([], c1);

% plotting
cfgp = [];
cfgp.layout    = 'EEG1010.lay';
cfgp.directionality = 'inflow';
cfgp.parameter = 'cohspctrm';
cfgp.refchannel = 'AFz';
ft_topoplotTFR(cfgp, c1);

% do JM's juicy transforms
tmp = ft_checkdata(freq, 'cmbrepresentation', 'fullfast');
c1             = ft_connectivityanalysis(cfgc, freq);
try
  ft_freqdescriptives([], tmp);
catch err
  fprintf('error: %s @ %s @ line %i\n ', err.message, err.stack(1).name, err.stack(1).line);
end

ft_topoplotTFR(cfgp, c1);
ft_connectivityplot([], c1);

cfga = [];
cfga.parameter = 'fourierspctrm';
try
  ft_appendfreq(cfga, freq, tmp);
catch err
  fprintf('error: %s @ %s @ line %i\n ', err.message, err.stack(1).name, err.stack(1).line);
end
vol = [];
vol.o = [0 0 4];
vol.r = 12;
vol.unit = 'cm';
vol.cond = 0.5;
vol.type = 'singlesphere';
cfgs.headmodel = vol;
cfgs.frequency = 10;
cfgs.elecfile = 'standard_1020.elc';
ft_sourceanalysis(cfgs, tmp);

tmp = ft_checkdata(tmp, 'cmbrepresentation', 'sparse');
c1             = ft_connectivityanalysis(cfgc, freq);
try
  ft_freqdescriptives([], tmp);
catch err
  fprintf('error: %s @ %s @ line %i\n ', err.message, err.stack(1).name, err.stack(1).line);
end
ft_topoplotTFR(cfgp, c1);
ft_connectivityplot([], c1);
try
  ft_appendfreq(cfga, freq, tmp);
catch err
  fprintf('error: %s @ %s @ line %i\n ', err.message, err.stack(1).name, err.stack(1).line);
end

vol = [];
vol.o = [0 0 4];
vol.r = 12;
vol.cond = 0.5;
vol.unit = 'cm';
vol.type = 'singlesphere';
cfgs.headmodel = vol;
cfgs.frequency = 10;
cfgs.elecfile = 'standard_1020.elc';
ft_sourceanalysis(cfgs, tmp);


tmp = ft_checkdata(tmp, 'cmbrepresentation', 'sparsewithpow');
c1             = ft_connectivityanalysis(cfgc, freq);
ft_freqdescriptives([], tmp);
ft_topoplotTFR(cfgp, c1);
ft_connectivityplot([], c1);
try
  ft_appendfreq(cfga, freq, tmp);
catch err
  fprintf('error: %s @ %s @ line %i\n ', err.message, err.stack(1).name, err.stack(1).line);
end

vol = [];
vol.o = [0 0 4];
vol.r = 12;
vol.cond = 0.5;
vol.unit = 'cm';
vol.type = 'singlesphere';
cfgs.headmodel = vol;
cfgs.elecfile = 'standard_1020.elc';
cfgs.frequency = 10;
ft_sourceanalysis(cfgs, tmp);

