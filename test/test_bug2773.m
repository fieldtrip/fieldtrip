function test_bug2773

% WALLTIME 00:20:00
% MEM 2gb

% TEST ft_dipolefitting ft_movieplotER ft_prepare_sourcemodel ft_prepare_layout


% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

orig = load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2773.mat'));
vol  = orig.cfg.vol;
elec = orig.cfg.elec;

cfg = [];
cfg.elec = elec;
cfg.channel =    {
  '1'
  '2'
  '3'
  '4'
  '5'
  '6'
  '7'
  '8'
  '9'
  '10'
  '11'
  '12'
  '13'
  '14'
  '15'
  '16'
  '17'
  '18'
  '19'
  '20'
  '21'
  '22'
  '23'
  '24'
  '25'
  '26'
  '27'
  '28'
  '29'
  '30'
  '31'
  '32'
  '33'
  '34'
  '35'
  '36'
  '37'
  '38'
  '39'
  '40'
  '41'
  '42'
  '43'
  '44'
  '45'
  '46'
  '47'
  '48'
  '49'
  '50'
  '51'
  '52'
  '53'
  '54'
  '55'
  '56'
  '57'
  '58'
  '59'
  '60'
  '61'
  };

layout = ft_prepare_layout(cfg); % it is not a very nice layout, but will do for the test script

figure; ft_plot_lay(layout);
figure; ft_plot_sens(elec);
figure; ft_plot_vol(vol);




