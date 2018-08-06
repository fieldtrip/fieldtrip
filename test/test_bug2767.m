function test_bug2767

% WALLTIME 00:10:00
% MEM 1gb


filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2767/01_ljh_firststd_meg_182.fif');

% this works, but returns hdr.nSamples=0
hdr = ft_read_header(filename);

try
  % this fails
  dat = ft_read_data(filename);
end

try
  % this fails
  evt = ft_read_event(filename);
end


% cfg = [];
% cfg.dataset = filename;
% raw = ft_preprocessing(cfg);

% cfg = [];
% erf = ft_timelockanalysis(cfg, raw);

