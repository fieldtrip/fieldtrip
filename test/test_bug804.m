function test_bug804

% WALLTIME 00:03:05

% TEST test_bug804
% TEST channeposition ft_datatype_sens yokogawa2grad ft_read_header

cd /home/common/matlab/fieldtrip/data/test/bug804

cfg = [];
cfg.dataset = 'amano.ave';
data = ft_preprocessing(cfg);

ft_plot_sens(data.grad);
