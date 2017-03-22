function test_bug2485

% WALLTIME 00:10:00
% MEM 500mb

% TEST ft_topoplotTFR

% Based on a script contributed by:
% PTB, T. Sander, 03072014
% tilmann.sander-thoemmes@ptb.de
% 
% Testscript for possible topoplot bug in fieldtrip-20140224
% filepath = 090512hulu.after1.con
% filename = YOUR PATH TO THE FILE

cfg             = [];
%cfg.dataset     = [filepath filename];
cfg.dataset     = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2485.con');
cfg.channel     = 'AG*';
cfg.trl= 1:100:900;
cfg.trl = [cfg.trl(:) cfg.trl(:)+99 zeros(length(cfg.trl), 1)];    
data            = ft_preprocessing(cfg); 

cfg = [];
cfg.channel={'all'};
cfg.keeptrials = 'yes';
trial_data = ft_timelockanalysis(cfg, data); 

cfg = [];
cfg.output ='powandcsd';
cfg.channel = [ trial_data.label(1:100)' trial_data.label(101:106)' ]';
for i = 1:100 
    cfg.channelcmb{i, 1}=trial_data.label{101};
    cfg.channelcmb{i, 2}=trial_data.label{i};
end
cfg.taper = 'hanning'; % 
cfg.method = 'mtmfft';
cfg.foilim = [0 70];
inp = ft_freqanalysis(cfg, trial_data);

cfg = [];
cfg.method = 'coh';
coh = ft_connectivityanalysis(cfg, inp);

cfg_topo         = [];
lay              = ft_prepare_layout(cfg_topo, data);
cfg_topo.layout  = lay;
ft_layoutplot(cfg_topo, data);       
cfg_topo.xlim= [0 5:10:75]; 
cfg_topo.refchannel = trial_data.label{101};
cfg_topo.parameter='cohspctrm';
cfg_topo.zlim=[-0.6 0.6]; 
ft_topoplotTFR(cfg_topo, coh);


