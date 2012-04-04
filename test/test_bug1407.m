function test_bug1407

% TEST test_bug1407
% TEST ft_read_header ft_read_data ft_read_event

% the following are from the fieldtrip/external/egi directory
% TEST read_mff_header read_mff_data read_mff_event mff_getEEGFilename mff_getSummaryInfo mff_getObject mff_micros2Sample

ft_hastoolbox('egi_mff', 1);
datadir = fileparts(which('read_mff_header'));
datadir = fullfile(datadir, 'test');
cd(datadir)

dataset = {
  'ATD256_3.ses.mff'
  'Long64ChannelWithEvents.mff'
  'NS1000sps.mff'
  'NS500Sine6Hz.mff'
  };

headerformat = 'egi_mff';
dataformat   = 'egi_mff';
eventformat  = 'egi_mff';

% headerformat = 'egi_mff_v1';
% dataformat   = 'egi_mff_v1';
% eventformat  = 'egi_mff_v1';

for i=1:length(dataset)
  hdr   = ft_read_header(dataset{i}, 'headerformat', headerformat);
  dat   = ft_read_data(dataset{i},   'headerformat', headerformat, 'dataformat', dataformat);
  event = ft_read_event(dataset{i},  'headerformat', headerformat, 'eventformat', eventformat);
  assert(size(dat,1)==length(hdr.label));
  assert(size(dat,2)==hdr.nSamples*hdr.nTrials);
end

