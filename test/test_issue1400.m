function test_issue1400

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_selectdata

ntrial = 2;

data = [];
data.label = {'1'};
for i=1:ntrial
  t = (1:1000)/1000 + i; % time vectors are unaligned to show the error
  s = sin(i*2*pi*t); % make a sine with 1, 2, 3, ... Hz
  data.time{i} = t;
  data.trial{i}(1,:) = s;  % insert it in the first channel
end

disp('Demonstrate averaging over time:')
disp('Averaging trials 1, 2')
cfg = [];
cfg.trials = [1,2];
cfg.avgoverrpt = 'yes';
cfg.keeprpt = 'no';
cfg.verbosity = 'off';
tmpdata = ft_selectdata(cfg, data);

% disp(['trialinfo(1, 1) = ', num2str(tmpdata.trialinfo(1,1))])
disp(['trial 1 time = ', regexprep(num2str(data.time{1}(1:5)),'\s+',', ')])
disp(['trial 2 time = ', regexprep(num2str(data.time{2}(1:5)),'\s+',', ')])
disp(['average time = ', regexprep(num2str(tmpdata.time{1}(1:5)),'\s+',', ')])

% add some trialinfo which is a mix of numbers and letters
disp('Demonstrate averaging over trialinfo:')
data.trialinfo = table([1;2],{'a'; 'b'});
disp(['trialinfo class is: ', class(data.trialinfo)])

% This would return an error in ft_selectdata, because it 
% averaging the table prior to discarding the averaged trialinfo
disp('Averaging trials 1, 2')
cfg = [];
cfg.trials = [1,2];
cfg.avgoverrpt = 'yes';
cfg.keeprpt = 'no';
cfg.verbosity = 'off';
tmpdata = ft_selectdata(cfg, data);
