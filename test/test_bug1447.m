function test_bug1447

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_multiplotER ft_singleplotER ft_plot_vector

% this contains the example data from Lilla, i.e. two ERPs and a layout
cd(dccnpath('/home/common/matlab/fieldtrip/data/test'))
load bug1447.mat

close all

% this can be a binary mask, or gradual between 0 and 1
% match.mask = abs(match.avg-nonmatch.avg)>2;

% use a gradial mask
% match.mask = ones(size(match.avg))/0.5;
% match.mask(1:30,801:1800) = 1;

% select some channels and some time for masking
match.mask = zeros(59, 2200);
match.mask(1:30,801:1800) = 1;
match.mask(15:end,401:600) = 1;

cfg = [];
cfg.layout = lay;
figure
ft_multiplotER(cfg, match, nonmatch);
title('no mask highlights');

% field in the first dataset to be used for marking significant data
cfg.maskparameter = 'mask';

% cfg.label = 'yes';
% cfg.axes  = 'x';
% cfg.box   = 'yes';
cfg.maskstyle     = 'box'; % style used for masking of data, 'box', 'thickness' or 'saturation' (default = 'box')
figure
ft_multiplotER(cfg, match, nonmatch);
title(sprintf('using %s highlights', cfg.maskstyle));

cfg.maskstyle     = 'thickness';
figure
ft_multiplotER(cfg, match, nonmatch);
title(sprintf('using %s highlights', cfg.maskstyle));

cfg.maskstyle     = 'saturation';
figure
ft_multiplotER(cfg, match, nonmatch);
title(sprintf('using %s highlights', cfg.maskstyle));

cfg.maskstyle     = 'difference';
figure
ft_multiplotER(cfg, match, nonmatch);
title(sprintf('using %s highlights', cfg.maskstyle));

