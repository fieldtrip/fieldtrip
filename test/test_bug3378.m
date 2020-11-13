function test_bug3378

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY

data = [];
data.powspctrm = randn(1,10,20);
data.time = 1:20;
data.freq = 1:10;
data.label = {'chan01'};
data.dimord = 'chan_freq_time';
data.mask = data.powspctrm>0;

figure('units','normalized','outerposition',[0 0 1 1]);

subplot(4,4,[5 6 9 10])
cfg=[];
cfg.colorbar = 'no';
ft_singleplotTFR(cfg,data);
title('Unmasked')

% case 1: logical mask, maskalpha < 1, opacity masking
cfg.maskparameter = 'mask';
cfg.maskalpha = 0.2;
subplot(4,4,3);
try
    ft_singleplotTFR(cfg, data); % this works as expected
catch ME
    disp(getReport(ME))
    text(0,.5, 'not supported')
    axis off
end
title('Logical, alpha~=1')

% case 2: numeric mask, maskalpha < 1, opacity masking
cfg.maskparameter = 'powspctrm';
cfg.maskalpha = 0.2;
subplot(4,4,4);
try
    ft_singleplotTFR(cfg, data); % this works as expected, warning is thrown
catch ME
    disp(getReport(ME))
    text(0,.5, 'not supported')
    axis off
end
title('Numeric, alpha~=1')

% case 3: logical mask, maskalpha = 1, opacity masking
cfg.maskparameter = 'mask';
cfg.maskalpha = 1;
subplot(4,4,7);
try
    ft_singleplotTFR(cfg, data); % this works as expected
catch ME
    disp(getReport(ME))
    text(0,.5, 'not supported')
    axis off
end
title('Logical, alpha=1')

% case 4: numeric mask, maskalpha = 1, opacity masking
cfg.maskparameter = 'powspctrm';
cfg.maskalpha = 1;
subplot(4,4,8);
try
    ft_singleplotTFR(cfg, data); % this works as expected, no warning with default maskalpha
catch ME
    disp(getReport(ME))
    text(0,.5, 'not supported')
    axis off
end
title('Numeric, alpha=1')

% case 5: logical mask, maskalpha < 1, outline masking
cfg.maskparameter = 'mask';
cfg.maskalpha = 0.2;
cfg.maskstyle = 'outline';
subplot(4,4,11);
try
    ft_singleplotTFR(cfg, data); % this works as expected - maskalpha is ignored
catch ME
    disp(getReport(ME))
    text(0,.5, 'not supported')
    axis off
end
title('Logical, alpha~=1')

% case 6: numeric mask, maskalpha < 1, outline masking
cfg.maskparameter = 'powspctrm';
cfg.maskalpha = 0.2;
subplot(4,4,12);
try
    ft_singleplotTFR(cfg, data); 
catch ME
    disp(getReport(ME))
    text(0,.5, 'not supported') % this works as expected
    axis off
end
title('Numeric, alpha~=1')

% case 7: logical mask, maskalpha = 1, outline masking
cfg.maskparameter = 'mask';
cfg.maskalpha = 1;
subplot(4,4,15);
try
    ft_singleplotTFR(cfg, data); % this works as expected - maskalpha is ignored
catch ME
    disp(getReport(ME))
    text(0,.5, 'not supported')
    axis off
end
title('Logical, alpha=1')

% case 8: numeric mask, maskalpha = 1, outline masking
cfg.maskparameter = 'powspctrm';
cfg.maskalpha = 1;
cfg.maskstyle = 'outline';
subplot(4,4,16);
try
    ft_singleplotTFR(cfg, data); 
catch ME
    disp(getReport(ME))
    text(0,.5, 'not supported') % this works as expected
    axis off
end
title('Numeric, alpha=1')

%% multiplot
load('easycapM25.mat')

data = [];
data.powspctrm = randn(23,10,20);
data.time = 1:20;
data.freq = 1:10;
data.label = lay.label;
data.dimord = 'chan_freq_time';
data.mask = data.powspctrm>0;

figure('units','normalized','outerposition',[0 0 1 1]);

subplot(4,4,[5 6 9 10])
cfg=[];
cfg.layout = lay;
cfg.colorbar = 'no';
ft_multiplotTFR(cfg,data);
title('Unmasked')

% case 1: logical mask, maskalpha < 1, opacity masking
cfg.maskparameter = 'mask';
cfg.maskalpha = 0.2;
subplot(4,4,3);
try
    ft_multiplotTFR(cfg, data); % this works as expected
catch ME
    disp(getReport(ME))
    text(0,.5, 'not supported')
    axis off
end
title('Logical, alpha~=1')

% case 2: numeric mask, maskalpha < 1, opacity masking
cfg.maskparameter = 'powspctrm';
cfg.maskalpha = 0.2;
subplot(4,4,4);
try
    ft_multiplotTFR(cfg, data); % this works as expected, warning is thrown
catch ME
    disp(getReport(ME))
    text(0,.5, 'not supported')
    axis off
end
title('Numeric, alpha~=1')

% case 3: logical mask, maskalpha = 1, opacity masking
cfg.maskparameter = 'mask';
cfg.maskalpha = 1;
subplot(4,4,7);
try
    ft_multiplotTFR(cfg, data); % this works as expected
catch ME
    disp(getReport(ME))
    text(0,.5, 'not supported')
    axis off
end
title('Logical, alpha=1')

% case 4: numeric mask, maskalpha = 1, opacity masking
cfg.maskparameter = 'powspctrm';
cfg.maskalpha = 1;
subplot(4,4,8);
try
    ft_multiplotTFR(cfg, data); % this works as expected, no warning with default maskalpha
catch ME
    disp(getReport(ME))
    text(0,.5, 'not supported')
    axis off
end
title('Numeric, alpha=1')

% case 5: logical mask, maskalpha < 1, outline masking
cfg.maskparameter = 'mask';
cfg.maskalpha = 0.2;
cfg.maskstyle = 'outline';
subplot(4,4,11);
try
    ft_multiplotTFR(cfg, data); % this works as expected - maskalpha is ignored
catch ME
    disp(getReport(ME))
    text(0,.5, 'not supported')
    axis off
end
title('Logical, alpha~=1')

% case 6: numeric mask, maskalpha < 1, outline masking
cfg.maskparameter = 'powspctrm';
cfg.maskalpha = 0.2;
subplot(4,4,12);
try
    ft_multiplotTFR(cfg, data); 
catch ME
    disp(getReport(ME))
    text(0,.5, 'not supported') % this works as expected
    axis off
end
title('Numeric, alpha~=1')

% case 7: logical mask, maskalpha = 1, outline masking
cfg.maskparameter = 'mask';
cfg.maskalpha = 1;
subplot(4,4,15);
try
    ft_multiplotTFR(cfg, data); % this works as expected - maskalpha is ignored
catch ME
    disp(getReport(ME))
    text(0,.5, 'not supported')
    axis off
end
title('Logical, alpha=1')

% case 8: numeric mask, maskalpha = 1, outline masking
cfg.maskparameter = 'powspctrm';
cfg.maskalpha = 1;
cfg.maskstyle = 'outline';
subplot(4,4,16);
try
    ft_multiplotTFR(cfg, data); 
catch ME
    disp(getReport(ME))
    text(0,.5, 'not supported') % this works as expected
    axis off
end
title('Numeric, alpha=1')
