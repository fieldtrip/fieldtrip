function failed_ft_datatype_source

% WALLTIME 24:00:00
% MEM 4gb

% TEST test_bug2185
% TEST ft_sourcegrandaverage ft_selectdata ft_selectdata_new ft_datatype_source ft_math

global ft_default
ft_default = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test the historical dirlists
% the style of this is also used in test_ft_analysisprotocol and test_ft_datatype


p = dccnpath('/home/common/matlab/fieldtrip/data/test/');
d = dir(fullfile(p, '20*'));
dirlist = {d.name};

% this cell-array will contain the filenames with data that were not converted
% according to the expectations
failure = {};

% these file contain data structures that have manually been identified as buggy,
% e.g. with inconsistent dimensions of the fields in source.trial.noisecsd
problematic = {
  dccnpath('/home/common/matlab/fieldtrip/data/test/20111231/source/meg/source_grid_mtmfft_fourier_trl_DICS_keepall_rawtrial_ctf275.mat')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20120630/source/meg/source_grid_mtmfft_fourier_trl_DICS_keepall_rawtrial_ctf275.mat')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20121231/source/meg/source_grid_mtmfft_fourier_trl_DICS_keepall_rawtrial_ctf275.mat')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20130630/source/meg/source_grid_mtmfft_fourier_trl_DICS_keepall_rawtrial_bti148.mat')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20130630/source/meg/source_grid_mtmfft_fourier_trl_DICS_keepall_rawtrial_ctf275.mat')
  };

for i=1:length(dirlist)
  % note that quite a few of the older versions do not have source data for testing
  d = dir(fullfile(p, dirlist{i}, 'source', 'meg', '*.mat'));
  filelist = {d.name};
  
  for j=1:length(filelist)
    
    filename = fullfile(p, dirlist{i}, 'source', 'meg', filelist{j});
    
    if any(strcmp(filename, problematic))
      fprintf('skipping buggy file %s\n', filename);
      continue
    end
    
    fprintf('file %d from %d, ', j, length(filelist));
    
    sourceold = loadvar(filename, 'source');
    sourcenew = ft_datatype_source(sourceold, 'dirlist', 'upcoming');
    
    fn = fieldnames(sourcenew);
    if all(cellfun(@isempty, strfind(fn, 'dim')))
      warning('this structure does not have a dimord field');
      failure{end+1} = filename;
    end
    
    clear sourceold sourcenew
  end % for j
end % for dirlist

fprintf('====================================================================================\n')
fprintf('the following historical source structures seem to cause condirlist problems\n')
fprintf('====================================================================================\n')
fprintf('%s\n', failure{:});
