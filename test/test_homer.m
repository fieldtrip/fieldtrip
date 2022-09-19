function test_homer

% WALLTIME 00:10:00
% MEM 2GB
% DEPENDENCY homer2fieldtrip fieldtrip2homer ft_write_data opto2homer homer2opto

% there are different ways that Homer data can be processed with FieldTrip
% - the way that is used in most FT tutorials is to read (original) Homer files using FT_PREPROCESSING
% - an alternative way is to read (original) Homer files using homer2fieldtrip
% - it is also possible to write FieldTrip data to a Homer file on disk

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/original/nirs/homer/S1001_run01.nirs');

data1 = homer2fieldtrip(filename);
data2 = homer2fieldtrip(fieldtrip2homer(data1));

p = tempdir;
f1 = fullfile(p, 'data1.nirs');
f2 = fullfile(p, 'data2.nirs');

ft_write_data(f1, data1.trial{1}, 'header', data1.hdr)
ft_write_data(f2, data2.trial{1}, 'header', data2.hdr)

% the files written by FieldTrip should still have the same content
data1a = homer2fieldtrip(f1);
data2a = homer2fieldtrip(f2);

% reading the files using the tutorial approach with FT_DEFINETRIAL and FT_PREPROCESSING
cfg = [];
cfg.dataset = filename;
data0b = ft_preprocessing(cfg);
cfg.dataset = f1;
data1b = ft_preprocessing(cfg);
cfg.dataset = f2;
data2b = ft_preprocessing(cfg);

% clean up
delete(f1)
delete(f2)

%%

% do not compare the hdr field, it is different in hdr.orig
fn = {'time', 'trial', 'label', 'opto'};

for i=1:numel(fn)
  field = fn{i};
  % after roundtrip fieldtrip->homer->fieldtrip
  assert(isalmostequal(data1.(field), data2. (field), 'abstol', 1e-8), sprintf('%s is different between data1 and data2', field));
  
  % after ft_write_data
  assert(isalmostequal(data1.(field), data1a.(field), 'abstol', 1e-8), sprintf('%s is different between data1 and data1a', field));
  assert(isalmostequal(data1.(field), data2a.(field), 'abstol', 1e-8), sprintf('%s is different between data1 and data2a', field));
  
  % after ft_preprocessing
  assert(isalmostequal(data1.(field), data0b.(field), 'abstol', 1e-8), sprintf('%s is different between data1 and data0b', field));
  assert(isalmostequal(data1.(field), data1b.(field), 'abstol', 1e-8), sprintf('%s is different between data1 and data1b', field));
  assert(isalmostequal(data1.(field), data2b.(field), 'abstol', 1e-8), sprintf('%s is different between data1 and data2b', field));
end
