function test_pull1856

% WALLTIME 00:20:00
% MEM 3gb
% DEPENDENCY ft_write_data getorthoviewpos

%%

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/original/nirs'));

%%

filename = {
  './artinis/Jorn/OM/1102 24 channel_2.oxy3'
  './artinis/Jorn/B24/TestData.oxy3'
  './artinis/Jorn/B24/B24_FT_S02.oxy3'
  './artinis/Jorn/B23/TestData.oxy3'
  './artinis/Helena/190529_walking.oxy3'
  % './artinis/Helena/sub-03_rec-01_nirs.oxy3' % this one does not have an optodetemplates.xml
  './artinis/Helena/190528_fingertap_L.oxy3'
  './artinis/Helena/190529_foottap_R.oxy3'
  };

%%

for i=1:numel(filename)
  cfg = [];
  cfg.dataset = filename{i};
  data_orig = ft_preprocessing(cfg);
  
  hdr = data_orig.hdr;
  dat = data_orig.trial{1};
  
  snirffile = [tempname '.snirf'];
  ft_write_data(snirffile, dat, 'header', hdr);
  
  cfg = [];
  cfg.dataset = snirffile;
  data_snirf = ft_preprocessing(cfg);
  
  delete(snirffile);
  
  % the channel names will be different, among others due to the rounding to the nominal wavelength
  % but at least the number of channels should be the same
  assert(length(data_orig.label)==length(data_snirf.label));
end
