function test_bug1785

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_read_sens read_asa_elc read_asa

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1785'));

filename = {
  'standard_1020.elc'
  'sphere_1020.elc'
  'waveguard64.elc'
  'EEG_Configuration_CP10040_Pilot13_Aug01_2012.txt'
  };

elec1 = ft_read_sens(filename{1}, 'fileformat', 'asa_elc');
elec2 = ft_read_sens(filename{2}, 'fileformat', 'asa_elc');
elec3 = ft_read_sens(filename{3}, 'fileformat', 'asa_elc');

% the fourth file turns out not to be an ASA file format, but a custom file format
% see http://bugzilla.fcdonders.nl/show_bug.cgi?id=1785#c4

% elec4 = ft_read_sens(filename{4}, 'fileformat', 'asa_elc');

