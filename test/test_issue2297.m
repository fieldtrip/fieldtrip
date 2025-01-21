filename = dccnpath('/project/3031000.02/test/original/meg/neuromag306/oddball_raw.fif');

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_trialfun_general
% DATA private

%%

% Maxshield data should be corrected using Maxfilter prior to importing in FieldTrip
% hence this should give an error

try
  ft_read_header(filename);
  ok = false;
catch
  ok = true;
end
assert(ok, 'this should have failed')

try
  ft_read_data(filename, 'begsample', 1, 'endsample', 1100);
  ok = false;
catch
  ok = true;
end
assert(ok, 'this should have failed')

try
  ft_read_event(filename);
  ok = false;
catch
  ok = true;
end
assert(ok, 'this should have failed')

%%

% this should work

ft_read_header(filename, 'checkmaxfilter', false);
ft_read_data(filename, 'begsample', 1, 'endsample', 1100, 'checkmaxfilter', false);
ft_read_event(filename, 'checkmaxfilter', false);

%%

% this should work

cfg = [];
cfg.dataset = filename;
cfg.trl = [1 1100 0];
cfg.checkmaxfilter = 'no';
ft_preprocessing(cfg);

%%

% this should work, but Pascal noticed that it failed
% due to the checkmaxfilter option not being passed down properly

cfg = [];
cfg.dataset = filename;
cfg.trialdef.eventtype = 'STI101';
cfg.trialdef.eventvalue = 64;
cfg.trialfun = 'ft_trialfun_general';
cfg.checkmaxfilter = 'no';
ft_definetrial(cfg);
