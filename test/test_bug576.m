function test_bug576

% TEST test_bug576
% TEST ft_checkdata ft_senstype

dataset122 = {
  '/home/common/matlab/fieldtrip/data/test/raw/meg/preproc_neuromag122.mat'
  '/home/common/matlab/fieldtrip/data/test/comp/meg/comp_neuromag122.mat'
  '/home/common/matlab/fieldtrip/data/test/freq/meg/freq_mtmfft_neuromag122.mat'
  '/home/common/matlab/fieldtrip/data/test/freq/meg/freq_mtmconvol_neuromag122.mat'
  '/home/common/matlab/fieldtrip/data/test/timelock/meg/timelock_neuromag122.mat'
  };

dataset306 = {
  '/home/common/matlab/fieldtrip/data/test/raw/meg/preproc_neuromag306.mat'
  '/home/common/matlab/fieldtrip/data/test/comp/meg/comp_neuromag306.mat'
  '/home/common/matlab/fieldtrip/data/test/freq/meg/freq_mtmfft_neuromag306.mat'
  '/home/common/matlab/fieldtrip/data/test/freq/meg/freq_mtmconvol_neuromag306.mat'
  '/home/common/matlab/fieldtrip/data/test/timelock/meg/timelock_neuromag306.mat'
  };


for i=1:length(dataset122)
  filename = dataset122{i};
  tmp = load(filename);
  fn = fieldnames(tmp);
  for j=1:length(fn)
    data = tmp.(fn{j});
    ft_checkdata(data, 'senstype', 'neuromag122');
  end % for j
end % for i


for i=1:length(dataset306)
  filename = dataset306{i};
  tmp = load(filename);
  fn = fieldnames(tmp);
  for j=1:length(fn)
    data = tmp.(fn{j});
    ft_checkdata(data, 'senstype', 'neuromag306');
  end % for j
end % for i

