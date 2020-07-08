function test_issue1167

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY edf2fieldtrip

filename = {'ma0844az_1-1+.edf';
            'test_generator.edf';
            'test_generator_2.edf'};
          
for k = 1:3
  data{k} = edf2fieldtrip(dccnpath(fullfile('/home/common/matlab/fieldtrip/data/test/original/eeg/edf', filename{k})));
end
