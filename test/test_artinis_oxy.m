function test_artinis_oxy

% WALLTIME 00:20:00
% MEM 4gb
% DEPENDENCY read_artinis_oxy3 read_artinis_oxy4

rootdir = dccnpath('/home/common/matlab/fieldtrip/data/test/original/nirs/artinis');

filename = {
  % oxy3
  'Jorn/OM/1102 24 channel_2.oxy3'
  'Jorn/B24/TestData.oxy3'
  'Jorn/B24/B24_FT_S02.oxy3'
  'Jorn/B23/TestData.oxy3'
  'Helena/190529_walking.oxy3'
  'Helena/190528_fingertap_L.oxy3'
  'Helena/190529_foottap_R.oxy3'
  % oxy4
  % 'Jorn/B23/TestData.oxy4'
  % 'Helena/190529_walking.oxy4'
  % 'Helena/190528_fingertap_L.oxy4'
  % 'Helena/190529_foottap_R.oxy4'
  };

% The oxy4 files result in the error
% [Artinis] Oxysoft not properly installed as a COM-Interface. Please see the appendix of the accompanying pdf-file!

%%
for i=1:numel(filename)
  
  file = fullfile(rootdir, filename{i});
  disp(file)
  
  [p, f, x] = fileparts(file);
  cd(p); % otherwise the GUI pops up, asking for the  optodetemplates.xml file
  
  hdr = ft_read_header(file);
  dat = ft_read_data(file);
  evt = ft_read_event(file);
end
