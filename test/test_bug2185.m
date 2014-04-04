function test_bug2185

% WALLTIME 00:20:00
% MEM 4gb

% TEST test_bug2185
% TEST ft_sourcegrandaverage ft_selectdata ft_selectdata_new ft_datatype_source ft_math

global ft_default
ft_default = [];

testsection1 = true;
testsection2 = true; % note that this takes a long time to test and will not throw an error in case something is wrong
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
  
  filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2185.mat');
  load(filename);
  
  % The content of the file is a cell-array that looks like this
  %
  % source_timelock_stim{1}
  % ans =
  %        time: [1x1500 double]
  %         pos: [8196x3 double]
  %      inside: [8196x1 double]
  %     outside: [1x0 double]
  %      method: 'average'
  %         avg: [1x1 struct]
  %         cfg: [1x1 struct]
  % source_timelock_stim{1}.avg
  % ans =
  %          mom: {8196x1 cell}
  %          pow: [8196x1500 double]
  %     noisecov: {8196x1 cell}
  
  cfg = [];
  cfg.parameter = 'pow';
  cfg.keepindividual = 'no';
  output = ft_sourcegrandaverage(cfg, source_timelock_stim{:})
  assert(isfield(output, 'pow'), 'missing output field');
  assert(isfield(output, 'time'), 'missing output field');
  
  cfg = [];
  cfg.parameter = 'pow';
  cfg.operation = 'multiply';
  cfg.value = 0;
  output = ft_math(cfg, output);
  assert(all(output.pow(:)==0), 'power should be zero');
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % intermezzo, do some interpolation and plotting
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  source2d = output;
  source2d.pow(:,:) = 0;
  source2d.pow(:,1) = source2d.pos(:,3); % replace the values by something easier to visualize
  
  % create a regular 3d-grid source structure that encompasses the cortical sheet
  minpos = floor(min(source2d.pos));
  maxpos = ceil(max(source2d.pos));
  xgrid = minpos(1):0.5:maxpos(1);
  ygrid = minpos(2):0.5:maxpos(2);
  zgrid = minpos(3):0.5:maxpos(3);
  dim = [length(xgrid) length(ygrid) length(zgrid)];
  [X, Y, Z] = ndgrid(xgrid, ygrid, zgrid);
  source3d.pos = [X(:) Y(:) Z(:)];
  source3d.dim = dim;
  % convert the regular 3d-grid source into a volume
  volume3d = ft_checkdata(source3d, 'datatype', 'volume')
  
  cfg = [];
  cfg.parameter = 'pow';
  source3d = ft_sourceinterpolate(cfg, source2d, volume3d)
  
  cfg = [];
  cfg.avgovertime = 'yes';
  cfg.keeptimedim = 'no';
  cfg.parameter = 'pow';
  source3davg = ft_selectdata(cfg, source3d)
  
  cfg = [];
  cfg.funparameter = 'pow';
  ft_sourceplot(cfg, source3davg);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % end of intermezzo
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  cfg = [];
  cfg.parameter = 'pow';
  cfg.keepindividual = 'yes';
  output = ft_sourcegrandaverage(cfg, source_timelock_stim{:})
  assert( isfield(output, 'pow'), 'missing output field');
  assert(isfield(output, 'time'), 'missing output field');
  
  cfg = [];
  cfg.parameter = 'pow';
  cfg.operation = 'multiply';
  cfg.value = 0;
  output = ft_math(cfg, output);
  assert(all(output.pow(:)==0), 'power should be zero');
  
  cfg = [];
  cfg.parameter = 'mom';
  cfg.keepindividual = 'no';
  output = ft_sourcegrandaverage(cfg, source_timelock_stim{:})
  assert( isfield(output, 'mom'), 'missing output field');
  assert( isfield(output, 'time'), 'missing output field');
  
  cfg = [];
  cfg.parameter = 'mom';
  cfg.operation = 'multiply';
  cfg.value = 0;
  output = ft_math(cfg, output);
  assert(all(output.mom{1}(:)==0), 'moment should be zero');
  
  cfg = [];
  cfg.parameter = 'mom';
  cfg.keepindividual = 'yes';
  output = ft_sourcegrandaverage(cfg, source_timelock_stim{:})
  assert( isfield(output, 'mom'), 'missing output field');
  assert( isfield(output, 'time'), 'missing output field');
  
  cfg = [];
  cfg.parameter = 'mom';
  cfg.operation = 'multiply';
  cfg.value = 0;
  output = ft_math(cfg, output);
  assert(all(output.mom{1}(:)==0), 'moment should be zero');
  
  cfg = [];
  cfg.parameter = 'noisecov';
  cfg.keepindividual = 'yes';
  output = ft_sourcegrandaverage(cfg, source_timelock_stim{:})
  
end % testsection3

