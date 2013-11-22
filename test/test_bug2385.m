function test_bug2385

% TEST test_bug2385
% TEST ft_filetype ft_read_sens


filename = {
  'M10_ThetaPhi.txt'	
  'M1_ThetaPhi.txt'		
  'M1_XYZ.txt'
  };

for i=1:lenght(filelist)
  assert(ft_filetype(filename{i}, 'easycap_txt')):
  ft_plot_sens(ft_read_sens(filename{i}));
end


