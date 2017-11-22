function test_bug1914

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_filetype ft_read_header ft_read_data ft_read_event

%%
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1914'));

dataset = {
  'conversion_testing_MEGEEG_raw.fif'
  'conversion_testing_MEGandperipheral_raw.fif'
  'conversion_testing_OnlyMEG_no_triggers_raw.fif'
  'conversion_testing_OnlyMEG_with_triggers_raw.fif'
  };

%%

numevent = [2 61 0 2421];

for i=1:length(dataset)
  filename = dataset{i};
  
  hdr   = ft_read_header(filename);
  dat   = ft_read_data(filename, 'begsample', 1, 'endsample', hdr.Fs);
  event = ft_read_event(filename);
  
  fprintf('============ %s ============\n', filename);
  fprintf('=======================================================================================\n')
  disp(hdr);
  disp(event);
  assert(length(event)==numevent(i));
  assert(ft_senstype(hdr, 'babysquid74'));  % the specific syste,
  assert(ft_senstype(hdr, 'babysquid'));    % more general, could also be artenis123
  
end
