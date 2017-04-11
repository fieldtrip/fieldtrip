function test_pull395

% WALLTIME 00:10:00
% MEM 2gb

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeg/bdf'));

%%

filename = {
  '050327BH_overCZnoAlpha.bdf'
  'Newtest17-2048.bdf'
  'Newtest17-256.bdf'
  'mbrain-train-smarting-mobile-eeg.bdf'
  };

for i=1:numel(filename)
  hdr = ft_read_header(filename{i});
  index = find(strcmpi(hdr.label, 'STATUS'));
  assert(hdr.orig.Cal(index)==1, 'calibration appears incorrect');
  
  sdat = ft_read_data(filename{i}, 'chanindx', index);
  assert(all((sdat-round(sdat))==0), 'non-integer values in STATUS channel');
  
  event = ft_read_event(filename{i});
  sel = strcmpi({event.type}, 'STATUS');
  smp = [event(sel).sample];
  val = [event(sel).value];
  figure; plot(smp, val, '.');
  assert(all((val-round(val))==0), 'non-integer values in STATUS events');
  assert(all(val>0), 'non-positive values in STATUS events');
  
end