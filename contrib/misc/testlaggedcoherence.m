datain=[];
datain.label={'Ch1','Ch2'};
datain.fsample=1000;
datain.time{1}=0:1/datain.fsample:1.832;
datain.trial{1}=randn(2,length(datain.time{1}));

cfg=[];
cfg.foi=[1 1.5 2 3 4 5 6 7 8 17 18 19 20 21];
cfg.loi=[1 2 3];
cfg.numcycles=3;

lcout=ft_laggedcoherence(cfg,datain);
