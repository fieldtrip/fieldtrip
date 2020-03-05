function test_pull1156

% WALLTIME 00:20:00
% MEM 3gb
% DEPENDENCY write_edf read_edf

%%

filename = {
  '/home/common/matlab/fieldtrip/data/test/original/eeg/edf/0601_s.edf'
  '/home/common/matlab/fieldtrip/data/test/original/eeg/edf/shhs1-200001.edf'
  '/home/common/matlab/fieldtrip/data/test/original/eeg/nicolet/Patient59_EEG-OPPTAKER-1_t1-EDF.edf'
  '/home/common/matlab/fieldtrip/data/test/original/eeg/nicolet/Patient59_EEG-OPPTAKER-1_t1-EDFPLUS.edf'
  '/home/common/matlab/fieldtrip/data/test/original/eeg/edf/RecordSession_2017.09.07_20.41.53.edf'
  };

headerfields = {
  'Fs'
  'nChans'
  'nTrials'
  'nSamples'
  'nSamplesPre'
  'label'
  };


%%

for i=1:numel(filename)
  hdr1 = ft_read_header(dccnpath(filename{i}));
  dat1 = ft_read_data(dccnpath(filename{i}));
  
  tempfile = [tempname '.edf'];
  ft_write_data(tempfile, dat1, 'header', hdr1);
  
  hdr2 = ft_read_header(tempfile);
  dat2 = ft_read_data(tempfile);
  
  delete(tempfile);
  
  if i~=5
    assert(isequal(keepfields(hdr1, headerfields), keepfields(hdr2, headerfields)));
    assert(corr(dat1(:), dat2(:))>0.99); % isalmostequal does not work well with all the different scales in the test files
  else
    % the 1st is in 1-sample blocks, whereas the 2nd is in 250 sample blocks
    % which causes the 2nd file to be slightly shorter
    hdr1.nSamples = 143000;
    hdr2.nSamples = 143000;
    dat1 = dat1(:,1:143000);
    dat2 = dat2(:,1:143000);
    assert(isequal(keepfields(hdr1, headerfields), keepfields(hdr2, headerfields)));
    assert(corr(dat1(:), dat2(:))>0.99); % isalmostequal does not work well with all the different scales in the test files
  end
end
