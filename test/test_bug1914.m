function test_bug1914

% MEM 1gb
% WALLTIME 00:03:17

% TEST test_bug1914
% TEST ft_filetype ft_read_header ft_read_data ft_read_event

cd(dccnfilename('/home/common/matlab/fieldtrip/data/test/bug1914'));

dataset = {
  'conversion_testing_MEGEEG_raw.fif'
  'conversion_testing_MEGandperipheral_raw.fif'
  'conversion_testing_OnlyMEG_no_triggers_raw.fif'
  'conversion_testing_OnlyMEG_with_triggers_raw.fif'
  };

for i=1:length(dataset)
  filename = dataset{i};
  
  hdr = ft_read_header(filename);
  dat = ft_read_data(filename, 'begsample', 1, 'endsample', hdr.Fs);
  event = ft_read_event(filename);
  
  fprintf('============ %s ============\n', filename);
  disp(hdr);
  disp(event);
  
end
