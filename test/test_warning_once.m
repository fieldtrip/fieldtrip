function test_warning_once

% MEM 1500mb
% WALLTIME 00:10:00

warning1 = 'hululu';
warning2 = 'aloah hey';

for i=1:2
  ft_warning('clear');
  ft_warning('once');
  ft_warning('backtrace', 'off');
  
  [output] = evalc('warning_caller(warning1, warning2)');
  w1size = strfind(output, warning1);
  w2size = strfind(output, warning2);
  if numel(w1size)<2 ||  numel(w2size)<1
    error('too few warnings thrown at iteration %d', i);
  elseif numel(w1size)>2 ||  numel(w2size)>1
    error('too many warnings thrown at iteration %d', i);
  end
end

% no clearing to verify whether these warnings stay!
[output] = evalc('warning_caller(warning1, warning2)');
w1size = strfind(output, warning1);
w2size = strfind(output, warning2);
if numel(w1size)<0 ||  numel(w2size)<0
  error('too few warnings thrown');
elseif numel(w1size)>0 ||  numel(w2size)>0
  error('too many warnings thrown');
end

% check some ft_ functions, therefore get dummy data
datainfo = ref_datasets;
dataset = datainfo(1);
load(fullfile(dataset.origdir, 'latest', 'raw', dataset.type, ['preproc_', dataset.datatype]));

cfg            = [];
cfg.method     = 'mtmconvol';
cfg.foi        = 1:.01:2;
cfg.taper      = 'hanning';
cfg.t_ftimwin  = ones(1, numel(cfg.foi)).*0.5;
cfg.toi        = (250:50:750)./1000;
cfg.polyremoval= 0;
% one warning should be thrown
ft_freqanalysis(cfg, data);
% again!
ft_freqanalysis(cfg, data);

end % function


function warning_caller(warning1, warning2)

for i=1:10
  ft_warning(warning1);
  ft_warning('FieldTrip:warning2', warning2);
end

% these warnings should be thrown now !again!
for i=1:10
  ft_warning(warning1);
  ft_warning('FieldTrip:warning2', warning2);
end

end % function

