function test_bug932

% MEM 2gb
% WALLTIME 00:20:00

% TEST ft_freqstatistics ft_selectdata ft_appendfreq ft_datatype_freq ft_datatype_sens

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug932'));

for tt=[4,6]
  load(strcat('LF_o_', num2str(tt)));
  load(strcat('LF_t_', num2str(tt)));
  sizeo = size(LF_o.powspctrm);
  sizet = size(LF_t.powspctrm);
  tatao= sizeo(1);
  tatat = sizet(1);
  
  cfg =[];
  cfg.channel = 'all';
  cfg.avgoverchan = 'yes';
  cfg.method = 'analytic';
  cfg.statistic = 'ft_statfun_indepsamplesT';
  cfg.design(2,1:tatao+tatat) = [1:tatao 1:tatat];
  cfg.design(1,1:tatao+tatat) = [ones(1,tatao) 2*ones(1,tatat)];
  cfg.ivar =1;
  eval(['t_t_o_subj' num2str(tt) ' = ft_freqstatistics(cfg,LF_o,LF_t);']);
end

