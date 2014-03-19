function test_bug2385

% TEST test_bug2385
% TEST ft_filetype ft_read_sens

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/original/electrodes/easycap'));

filename = {
  'M10_ThetaPhi.txt'
  'M1_ThetaPhi.txt'
  'M1_XYZ.txt'
  };

for i=1:length(filename)
  assert(ft_filetype(filename{i}, 'easycap_txt'));
  elec = ft_read_sens(filename{i});
  figure
  ft_plot_sens(elec);
end


