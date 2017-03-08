function test_bug3078

% MEM 500mb
% WALLTIME 00:15:00

% TEST ft_sourceanalysis

% there's something obscure going on with the channel order in the spatial filters,
% supposedly in an interaction with precomputed leadfields and/or freq data with powandcsd

% create some data
data = [];
data.label = {'d';'c';'b';'a'}; 
% do 4 channels: when playing with 3 channels and exploring the effect of 
% taking a subset of channels in the cfg, I noticed strange behaviour when
% converting sparse->full, or some interaction with ft_selectdata (I
% suspect the haschancmb situation not being handled well with
% 'sparsewithpow' input, if the number of combinations is equal to the
% number of labels. this is very specifically the case if the input data
% has a complete sparsewithpow representation with 3 channels (yielding 3
% labelcmbs), as well as three labels. FIXME: THIS SHOULD BE FOLLOWED UP
% ELSEWHERE
for k = 1:5
	data.trial{k} = randn(4,1000);
	data.time{k}  = (0:999)./1000;
end
data.grad.coilpos = [0 0 8;0 -8 0;0 8 0;8 0 0];
data.grad.coilori = [0 0 1;0 -1 0;0 1 0;1 0 0];
data.grad.tra     = eye(4);
data.grad.label   = {'a';'b';'c';'d'};
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
cfg.dics.keepfilter    = 'yes';
cfg.dics.realfilter    = 'yes';
source1 = ft_sourceanalysis(cfg, freq1);
source2 = ft_sourceanalysis(cfg, freq2);

firstinside         = find(source1.inside,1,'first');
channelorderflipped = (isalmostequal(fliplr(source1.avg.filter{firstinside}),source2.avg.filter{firstinside},'reltol',1e-4));
if channelorderflipped,
	fprintf('Currently the order of the channels in the spatial filter is flipped, depending on whether ''fourier'' or ''powandcsd'' in input\n');
  % this is problematic if true, because at present there's no information
  % with respect to the order of the channels in the output...
end
 
% if channelorderflipped==0, there's no issue anymore with the freq1/freq2
% channel order, let's proceed with checking the effect of channel order in the
% leadfields

cfg = [];
cfg.method = 'dics';
cfg.frequency = 20;
cfg.headmodel = headmodel;
cfg.grid      = sourcemodel_lf1;
cfg.dics.keepleadfield = 'yes';
cfg.dics.keepfilter    = 'yes';
cfg.dics.realfilter    = 'yes';
source1 = ft_sourceanalysis(cfg, freq1);
cfg.grid      = sourcemodel_lf2;
source2 = ft_sourceanalysis(cfg, freq1);

firstinside         = find(source1.inside,1,'first');
channelordernotflipped = (isalmostequal(source1.avg.filter{firstinside},source2.avg.filter{firstinside},'reltol',1e-4));
if channelordernotflipped,
	fprintf('Currently the order of the channels in the spatial filter is identical, despite the channel ordering in the input leadfields\n');
  % this is good, because it suggests that the order of labels in the data determines the channel order...
	% what's not good, is that there's now way to know the order for sure,
	% and it will be wrong if it's assumed to be the order in the leadfield
end


% now see what happens with different settings of keepleadfield
cfg = [];
cfg.method = 'dics';
cfg.frequency = 20;
cfg.headmodel = headmodel;
cfg.grid      = sourcemodel;
cfg.dics.keepfilter    = 'yes';
cfg.dics.realfilter    = 'yes';
source1 = ft_sourceanalysis(cfg, freq1);
cfg.keepleadfield = 'yes';
source2 = ft_sourceanalysis(cfg, freq1);
cfg.grid      = sourcemodel_lf1;
source3 = ft_sourceanalysis(cfg, freq1);

