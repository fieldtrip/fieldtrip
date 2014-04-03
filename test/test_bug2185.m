% function test_bug2185

% WALLTIME 00:20:00
% MEM 4gb

% TEST test_bug2185
% TEST ft_sourcegrandaverage ft_selectdata ft_selectdata_new ft_datatype_source

global ft_default
ft_default = [];

testsection1 = false;
testsection2 = false;
testsection3 = true;

if testsection1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% use some hand-crafted data, version 1
  
  source = [];
  source.dim = [10 11 12];
  source.transform = eye(4);
  source.avg.pow = rand(10*11*12,1);
  source.inside = 1:660;
  source.outside = 661:1320;
  
  ft_checkdata(source, 'datatype', 'source')
  
  cfg = [];
  cfg.parameter = 'pow';
  cfg.keepindividual = 'no';
  grandavg = ft_sourcegrandaverage(cfg, source, source)
  
  cfg.keepindividual = 'yes';
  grandavg = ft_sourcegrandaverage(cfg, source, source)
  
  %% version 2
  
  source = [];
  source.transform = eye(4);
  source.pos = rand(10,3);
  source.pow = rand(10,1);
  source.inside = 1:660;
  source.outside = 661:1320;
  
  
  ft_checkdata(source, 'datatype', 'source')
  
  cfg = [];
  cfg.parameter = 'pow';
  cfg.keepindividual = 'no';
  grandavg = ft_sourcegrandaverage(cfg, source, source)
  
  cfg.keepindividual = 'yes';
  grandavg = ft_sourcegrandaverage(cfg, source, source)
  
  %% version 3
  
  source = [];
  source.dim = [10 11 12];
  source.transform = eye(4);
  source.pow = rand(10*11*12,1);
  source.inside = 1:660;
  source.outside = 661:1320;
  
  ft_checkdata(source, 'datatype', 'source')
  
  cfg = [];
  cfg.parameter = 'pow';
  cfg.keepindividual = 'no';
  grandavg = ft_sourcegrandaverage(cfg, source, source)
  
  cfg.keepindividual = 'yes';
  grandavg = ft_sourcegrandaverage(cfg, source, source)
  
end % testsection1

if testsection2
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% test the historical versions
  clear all
  
  p = dccnpath('/home/common/matlab/fieldtrip/data/test/');
  d = dir(fullfile(p, '20*'));
  version = {d.name};
  
  failure = {};
  
  for i=1:length(version)
    % note that quite a few of the older versinos don't have source data for testing
    d = dir(fullfile(p, version{i}, 'source', 'meg', '*.mat'));
    filelist = {d.name};
    for j=1:length(filelist)
      
      filename = fullfile(p, version{i}, 'source', 'meg', filelist{j});
      varname = 'source';
      
      % these file contain data structures that have manually been identified as buggy,
      % e.g. with inconsistent dimensions of the fields in source.trial.noisecsd
      problematic = {
        '/home/common/matlab/fieldtrip/data/test/20111231/source/meg/source_grid_mtmfft_fourier_trl_DICS_keepall_rawtrial_ctf275.mat'
        '/home/common/matlab/fieldtrip/data/test/20120630/source/meg/source_grid_mtmfft_fourier_trl_DICS_keepall_rawtrial_ctf275.mat'
        '/home/common/matlab/fieldtrip/data/test/20121231/source/meg/source_grid_mtmfft_fourier_trl_DICS_keepall_rawtrial_ctf275.mat'
        '/home/common/matlab/fieldtrip/data/test/20130630/source/meg/source_grid_mtmfft_fourier_trl_DICS_keepall_rawtrial_bti148.mat'
        '/home/common/matlab/fieldtrip/data/test/20130630/source/meg/source_grid_mtmfft_fourier_trl_DICS_keepall_rawtrial_ctf275.mat'
        };
      
      if any(strcmp(filename, problematic))
        fprintf('skipping buggy file %s\n', filename);
        continue
      end
      
      fprintf('processing %s\n', filename);
      
      %%%% this section is from private/loadvar
      var = whos('-file', filename);
      if length(var)==1
        filecontent = load(filename); % read the one variable in the file, regardless of how it is called
        value       = filecontent.(var.name);
        clear filecontent
      else
        filecontent = load(filename, varname);
        value       = filecontent.(varname);  % read the variable named according to the input specification
        clear filecontent
      end
      %%%% end of section from private/loadvar
      sourceold = value;
      clear value
      
      sourcenew = ft_datatype_source(sourceold, 'version', 'upcoming');
      
      fn = fieldnames(sourcenew);
      if all(cellfun(@isempty, strfind(fn, 'dim')))
        warning('this structure does not have a dimord field');
        failure{end+1} = filename;
      end
      
      clear sourceold sourcenew
    end % for j
  end % for version
  
  fprintf('====================================================================================\n')
  fprintf('the following historical source structures seem to cause conversion problems\n')
  fprintf('====================================================================================\n')
  fprintf('%s\n', failure{:});
  
end % testsection2

if testsection3
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% use the end-user provided test data
  clear all
  
  filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2185.mat');
  load(filename);
  
  cfg = [];
  cfg.parameter = 'pow';
  cfg.keepindividual = 'no';
  output = ft_sourcegrandaverage(cfg, source_timelock_stim{:});
  
  assert( isfield(output, 'pow'), 'missing output field');
  assert(~isfield(output, 'time'), 'the output should not have time, only pow');
  
  cfg = [];
  cfg.parameter = 'pow';
  cfg.keepindividual = 'yes';
  output = ft_sourcegrandaverage(cfg, source_timelock_stim{:});
  
  assert( isfield(output, 'pow'), 'missing output field');
  assert(~isfield(output, 'time'), 'the output should not have time, only pow');
  
  cfg = [];
  cfg.parameter = 'mom';
  cfg.keepindividual = 'no';
  output = ft_sourcegrandaverage(cfg, source_timelock_stim{:});
  
  assert( isfield(output, 'mom'), 'missing output field');
  assert( isfield(output, 'time'), 'missing output field');
  
  cfg = [];
  cfg.parameter = 'mom';
  cfg.keepindividual = 'yes';
  output = ft_sourcegrandaverage(cfg, source_timelock_stim{:});
  
  assert( isfield(output, 'mom'), 'missing output field');
  assert( isfield(output, 'time'), 'missing output field');

    cfg = [];
  cfg.parameter = 'noisecov';
  cfg.keepindividual = 'yes';
  output = ft_sourcegrandaverage(cfg, source_timelock_stim{:});

  
end % testsection3

