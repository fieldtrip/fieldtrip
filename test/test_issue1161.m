function test_issue1161

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY snirf

%%

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/original/nirs/snirf'));

filename = {
  'SNIRF_example.snirf'
  'Simple_Probe.snirf'
  'neuro_run01.snirf'
  };

%%

for i=1:numel(filename)
  hdr = ft_read_header(filename{i});
  dat = ft_read_data(filename{i});
  evt = ft_read_event(filename{i});
end
