function test_ft_dipolefitting

% WALLTIME 00:20:00
% MEM 2gb
% DEPENDENCY ft_dipolefitting ft_headmodel_concentricspheres
% DATATYPE comp timelock freq

fs = 500;
nchan = 32;
start_time = -1; %s
end_time = 2.5; %s
nsamples = (end_time - start_time) * fs + 1;

data = [];
data.time{1} = linspace(start_time, end_time, nsamples);
data.trial{1} = randn(nchan,nsamples);
data.sampleinfo = [1 nsamples];
data.label = cellstr(num2str((1:nchan).'));
data.sens.label = data.label;
data.sens.chanpos = randn(nchan,3);
data.sens.elecpos = data.sens.chanpos;
data.sens.tra = eye(nchan);
data.sens.unit = 'cm';

geometry = [];
geometry.pos = data.sens.chanpos;
geometry.unit = data.sens.unit;
cfg = [];
cfg.conductivity = [0.33 1.00 0.042 0.33];
headmodel = ft_headmodel_concentricspheres(geometry, 'conductivity', cfg.conductivity);

% test for comp type
cfg = [];
comp = ft_componentanalysis(cfg,data);
ft_checkdata(comp, 'datatype', 'comp');
cfg = [];
cfg.xgrid = 'auto';
cfg.ygrid = 'auto';
cfg.zgrid = 'auto';
cfg.resolution = 1; %cm
cfg.component = 1:5;
cfg.headmodel = headmodel;
cfg.elec = data.sens;
sourceout = ft_dipolefitting(cfg, comp);

% test for timelock type
cfg = [];
timelock = ft_timelockanalysis(cfg,data);
ft_checkdata(timelock, 'datatype', 'timelock');
cfg = [];
cfg.xgrid = 'auto';
cfg.ygrid = 'auto';
cfg.zgrid = 'auto';
cfg.resolution = 1; %cm
cfg.latency = [0 1];
cfg.model = 'regional';
cfg.headmodel = headmodel;
cfg.elec = data.sens;
sourceout = ft_dipolefitting(cfg, timelock);

% test for freq type
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.tapsmofrq = 2;
cfg.output = 'fourier';
freq = ft_freqanalysis(cfg,data);
ft_checkdata(freq, 'datatype', 'freq');
cfg = [];
cfg.xgrid = 'auto';
cfg.ygrid = 'auto';
cfg.zgrid = 'auto';
cfg.resolution = 1; %cm
cfg.frequency = 10; %Hz
cfg.headmodel = headmodel;
cfg.elec = data.sens;
sourceout = ft_dipolefitting(cfg, freq);
