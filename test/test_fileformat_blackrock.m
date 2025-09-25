function test_fileformat_blackrock

% MEM 4gb
% WALLTIME 00:20:00
% DEPENDENCY ft_read_data ft_read_header ft_read_event
% DATA private

% check the data in the 'lfp' folder
datadir = dccnpath('/project/3031000.02/test/original/lfp/blackrock');

datasets = {'FileSpec2.3001.ccf'
    'FileSpec2.3001.nev'
    'FileSpec2.3001.ns5'
    'FirstTest001.ccf'
    'FirstTest001.nev'
    'FirstTest001.ns3'
    'FirstTest001.ns6'
    'NS4Test001.ccf'
    'NS4Test001.nev'
    'NS4Test001.ns4'
    'NS4Test001.ns6'
    'SecondTest002.ccf'
    'SecondTest002.nev'
    'SecondTest002.ns3'
    'SecondTest002.ns6'};



ok = true(numel(datasets),1);

for k = 1:numel(datasets)
  cfg = [];
  cfg.dataset = fullfile(datadir, datasets{k});
  try
    if endsWith(datasets{k}, 'nev')
      data = ft_read_spike(cfg.dataset);
    else
      data = ft_preprocessing(cfg);
    end
  catch ME
    err{k,1} = ME;
    
    ok(k,1) = false;
  end
end

fprintf('\n');
% it looks like files with the extension ccf are not supported
for k = find(~ok)'
  fprintf('failed to read file %s\n', datasets{k});
  fprintf('%s\n', err{k}.message);
  fprintf('%s\n', err{k}.identifier);
  fprintf('\n');
end


% also check the test data in the 'eeg' folder
datadir = dccnpath('/project/3031000.02/test/original/eeg/blackrock');

datasets = {'LT1D0.000F0000.mat'
  'test_data_blackrock.ns2'};

cfg.dataset = fullfile(datadir, datasets{1});
data = ft_preprocessing(cfg);
cfg.dataset = fullfile(datadir, datasets{2});
data = ft_preprocessing(cfg);


