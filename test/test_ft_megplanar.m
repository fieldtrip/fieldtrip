function test_ft_megplanar

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_megplanar ft_prepare_neighbours ft_topoplotER

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

datainfo = ref_datasets;
sel      = match_str({datainfo.datatype}',{'bti148' 'bti248' 'ctf151' 'ctf275' 'itab153' 'yokogawa160'}');
datainfo = datainfo(sel);

%% check whether all runs through smoothly
for k = 1:numel(datainfo)
  fname = fullfile(datainfo(k).origdir, 'latest/raw',datainfo(k).type,['preproc_' datainfo(k).datatype]);
  load(fname);

  cfg = [];
  cfg.demean = 'yes';
  cfg.feedback = 'textbar';
  data = ft_preprocessing(cfg, data);
  
  cfg = [];
  cfg.method = 'triangulation';
  cfg.channel = 'MEG';
  neighbours = ft_prepare_neighbours(cfg, data);

  cfg = [];
  cfg.neighbours   = neighbours;
  cfg.planarmethod = 'sincos';
  data2 = ft_megplanar(cfg, data);
  
  cfg.planarmethod = 'orig';
  data3 = ft_megplanar(cfg, data);
  
  cfg.planarmethod = 'fitplane';
  data4 = ft_megplanar(cfg, data);
  
  vol      = [];
  vol.o    = [0 0 4];
  vol.r    = 12;
  vol.unit = 'cm';
  
  % ensure units in the gradiometer array and volume conductor to be equal
  data.grad = ft_convert_units(data.grad);
  vol       = ft_convert_units(vol, data.grad.unit);
  cfg.planarmethod = 'sourceproject';
  cfg.vol = vol;
  data5 = ft_megplanar(cfg, data);
  
%   % subtract mean
%   cfg = [];
%   cfg.demean   = 'yes';
%   cfg.feedback = 'textbar';
%   data2 = ft_preprocessing(cfg, data2);
%   data3 = ft_preprocessing(cfg, data3);
%   data4 = ft_preprocessing(cfg, data4);
%   data5 = ft_preprocessing(cfg, data5);
%   
%   cfg = [];
%   data2 = ft_combineplanar(cfg, data2);
%   data3 = ft_combineplanar(cfg, data3);
%   data4 = ft_combineplanar(cfg, data4);
%   data5 = ft_combineplanar(cfg, data5);
%   
%   % subtract mean
%   cfg = [];
%   cfg.demean   = 'yes';
%   cfg.feedback = 'textbar';
%   data2 = ft_preprocessing(cfg, data2);
%   data3 = ft_preprocessing(cfg, data3);
%   data4 = ft_preprocessing(cfg, data4);
%   data5 = ft_preprocessing(cfg, data5);
%   
%   % compute power of channels
%   pow2 = sum(cat(2,data2.trial{:}).^2,2);
%   pow3 = sum(cat(2,data3.trial{:}).^2,2);
%   pow4 = sum(cat(2,data4.trial{:}).^2,2);
%   pow5 = sum(cat(2,data5.trial{:}).^2,2);
%  
%   % compute covariance between same channels
%   cov23 = sum(cat(2,data2.trial{:}).*cat(2,data3.trial{:}),2);
%   %cov24 = sum(cat(2,data2.trial{:}).*cat(2,data4.trial{:}),2);
%   cov25 = sum(cat(2,data2.trial{:}).*cat(2,data5.trial{:}),2);
%   %cov34 = sum(cat(2,data3.trial{:}).*cat(2,data4.trial{:}),2);
%   cov35 = sum(cat(2,data3.trial{:}).*cat(2,data5.trial{:}),2);
%   %cov45 = sum(cat(2,data4.trial{:}).*cat(2,data5.trial{:}),2);
%   
%   % compute correlation
%   corr23 = cov23./sqrt(pow2.*pow3);
%   %corr24 = cov24./sqrt(pow2.*pow4);
%   corr25 = cov25./sqrt(pow2.*pow5);
%   %corr34 = cov34./sqrt(pow3.*pow4);
%   corr35 = cov35./sqrt(pow3.*pow5);
%   %corr45 = cov45./sqrt(pow4.*pow5);
  
end

%% check for qualitative changes in transformation
% load exampledata from ER tutorial, simplified

cd(dccnpath('/home/common/matlab/fieldtrip/data/test'))
load avgFIC

% compute neighbours
cfg = [];
cfg.method =  'distance';
cfg.neighbours = ft_prepare_neighbours(cfg, avgFIC);

% megplanar
cfg.planarmethod = 'sincos';
avgFICplanar = ft_megplanar(cfg, avgFIC);

% timelock and combine planar
cfg = [];
avgFICplanar     = ft_timelockanalysis(cfg, avgFICplanar);
avgFICplanarComb = ft_combineplanar(cfg,avgFICplanar);


% check whether maxima and minima on sensor level are the same as in summer
% 2011
it = 5; % #iterations
i = 0;
while it < i
    [a, b] = nanmax(nanmax(avgFICplanarComb.avg, [], 2));
    if (b~=29) || b~= 53 || b~=4 || b~=28 || b~=3
        error('the global maxima has moved location');
    else
        avgFICplanarComb.avg(b, :) = nan;
    end
    i = t+1;
end
fprintf('global maxima seems all fine\n');

i = 0;
while it < i
    [a, b] = nanmin(nanmin(avgFICplanarComb.avg, [], 2));
    if b~=21|| b~=33 || b~=133 || b~=92 || b~=144
        error('the global minima has moved location');
    else
        avgFICplanarComb.avg(b, :) = nan;
    end
    i = t+1;
end
fprintf('global minima seems all fine\n');


% plotting
cfg = [];
clf
subplot(121)
cfg.xlim = [0.3 0.5];
cfg.layout = 'CTF151.lay';
cfg.zlim = 'maxmin';
cfg.colorbar = 'yes';
ft_topoplotER(cfg,avgFIC)
colorbar
subplot(122)
cfg.zlim = 'maxabs';
% cfg.layout = 'CTF151.lay';
ft_topoplotER(cfg,avgFICplanarComb)

fprintf('compare with http://www.fieldtriptoolbox.org/tutorial/eventrelatedaveraging#plot_the_results_planar_gradients\n');
