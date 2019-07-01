function test_pull1156

% WALLTIME 00:20:00
% MEM 2gb
% DEPENDENCY write_edf read_edf

%%

filename = {
  '/home/common/matlab/fieldtrip/data/test/original/eeg/edf/0601_s.edf'
  '/home/common/matlab/fieldtrip/data/test/original/eeg/edf/shhs1-200001.edf'
  '/home/common/matlab/fieldtrip/data/test/original/eeg/edf/RecordSession_2017.09.07_20.41.53.edf'
  '/home/common/matlab/fieldtrip/data/test/original/eeg/nicolet/Patient59_EEG-OPPTAKER-1_t1-EDF.edf'
  '/home/common/matlab/fieldtrip/data/test/original/eeg/nicolet/Patient59_EEG-OPPTAKER-1_t1-EDFPLUS.edf'
  };

headerfields = {
  'Fs'
  'nChans'
  'label'
  'nSamples'
  'nSamplesPre'
  'nTrials'
  };

%%

for i=1:numel(filename)
  hdr = ft_read_header(dccnpath(filename{i}));
  dat = ft_read_data(dccnpath(filename{i}));
  
  tempfile = [tempname '.edf'];
  ft_write_data(tempfile, dat, 'header', hdr);
  
  hdr1 = ft_read_header(tempfile);
  dat1 = ft_read_data(tempfile);
  
  assert(isequal(keepfields(hdr, headerfields), keepfields(hdr1, headerfields)));
  assert(corr(dat(:), dat1(:))>0.99); % isalmostequal does not work well with all the different scales in the test files
end
