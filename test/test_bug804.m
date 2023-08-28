function test_bug804

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY channeposition ft_datatype_sens yokogawa2grad ft_read_header
% DATA private

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug804'));

cfg = [];
cfg.dataset = 'amano.ave';
data = ft_preprocessing(cfg);

ft_plot_sens(data.grad);
