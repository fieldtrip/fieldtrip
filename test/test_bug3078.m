function test_bug3078

% MEM=500mb
% WALLTIME=00:15:00

% TEST test_bug3078
% TEST ft_sourceanalysis

% there's something obscure going on with the channel order in the spatial filters,
% supposedly in an interaction with precomputed leadfields and/or freq data with powandcsd

% create some data
data = [];
data.label = {'c';'b';'a'};
for k = 1:5
	data.trial{k} = randn(3,1000);
	data.time{k}  = (0:999)./1000;
end
data.grad.coilpos = [0 0 8;0 -8 0;0 8 0];
data.grad.coilori = [0 0 1;0 -1 0;0 1 0];
data.grad.tra     = eye(3);
data.grad.label   = {'a';'b';'c'};
data.grad.chanpos = data.grad.coilpos;
data.grad.chanori = data.grad.coilori;
data.grad = ft_datatype_sens(data.grad);

cfg = [];
cfg.method = 'mtmfft';
cfg.foilim = [20 20];
cfg.output = 'powandcsd';
cfg.channelcmb = {'all' 'all'};
cfg.tapsmofrq = 2;
freq1 = ft_freqanalysis(cfg, data);

cfg   = rmfield(cfg,'channelcmb');
cfg.output = 'fourier';
freq2 = ft_freqanalysis(cfg, data);

headmodel.o = [0 0 0];
headmodel.r = 8;
headmodel.type = ft_voltype(headmodel);
headmodel   = ft_datatype_headmodel(headmodel);

cfg = [];
cfg.grid.resolution = 2;
cfg.headmodel = headmodel;
cfg.grad      = data.grad;
sourcemodel   = ft_prepare_sourcemodel(cfg);

% compute the leadfields with a grad in the cfg, returns order as in the
% grad.label
cfg           = [];
cfg.grid      = sourcemodel;
cfg.headmodel = headmodel;
cfg.grad      = data.grad;
sourcemodel_lf1 = ft_prepare_leadfield(cfg);

% compute the leadfields with a data argument, returns order as in
% data.label
cfg             = rmfield(cfg, 'grad');
sourcemodel_lf2 = ft_prepare_leadfield(cfg, data);

cfg = [];
cfg.method = 'dics';
cfg.frequency = 20;
cfg.headmodel = headmodel;
cfg.grid      = sourcemodel;
cfg.dics.keepleadfield = 'yes';
cfg.dics.keepfilter    = 'yes';
cfg.dics.realfilter    = 'yes';
source1 = ft_sourceanalysis(cfg, freq1);
source2 = ft_sourceanalysis(cfg, freq2);

firstinside = find(source1.inside,1,'first');
channelorderflipped = (identical(fliplr(source1.avg.filter{firstinside}),source2.avg.filter{firstinside},'reltol',1e-4))
if channelorderflipped,
	fprintf('Currently the order of the channels in the spatial filter is flipped, depending on whether ''fourier'' or ''powandcsd'' in input\n');
  % this is problematic if true, because at present there's no information
  % with respect to the order of the channels in the output...
end
