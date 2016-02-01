function test_bug1248

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug1248 ft_preprocessing preproc

for samp=[25 98 999]
  for fs=[235 943]
    data = [];
    data.hdr.Fs     = fs;
    data.hdr.label  = {'fakechan1','fakechan2','fakechan3'};
    data.fsample    = fs;
    data.label      = {'fakechan1','fakechan2','fakechan3'};
    data.trial      = {rand(3,samp),rand(3,samp+1),rand(3,samp+2)};
    data.time       = {1/fs:1/fs:samp/fs,1/fs:1/fs:(samp+1)/fs,1/fs:1/fs:(samp+2)/fs};
    
    cfg = [];
    cfg.hpfreq     = 2;
    cfg.hpfiltord  = [];
    cfg.hpfilter   = 'yes';
    cfg.hpfilttype = 'fir';
    data = ft_preprocessing(cfg,data);
    
    cfg = [];
    cfg.hpfreq     = 2;
    cfg.hpfiltord  = [];
    cfg.hpfilter   = 'yes';
    cfg.hpfilttype = 'firls';
    data = ft_preprocessing(cfg,data);
    
    % no problem originally here
    cfg = [];
    cfg.lpfreq     = 20;
    cfg.lpfiltord  = [];
    cfg.lpfilter   = 'yes';
    cfg.lpfilttype = 'fir';
    data = ft_preprocessing(cfg,data);
    
    cfg = [];
    cfg.bpfreq     = [2 20];
    cfg.bpfiltord  = [];
    cfg.bpfilter   = 'yes';
    cfg.bpfilttype = 'fir';
    data = ft_preprocessing(cfg,data);
    
    cfg = [];
    cfg.bsfreq     = [2 20];
    cfg.bsfiltord  = [];
    cfg.bsfilter   = 'yes';
    cfg.bsfilttype = 'fir';
    data = ft_preprocessing(cfg,data);
    
    cfg = [];
    cfg.bsfreq     = [2 20];
    cfg.bsfiltord  = [];
    cfg.bsfilter   = 'yes';
    cfg.bsfilttype = 'firls';
    data = ft_preprocessing(cfg,data);
    
  end
end
