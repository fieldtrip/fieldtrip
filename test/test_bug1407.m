function test_bug1407

% TEST test_bug1407
% TEST ft_read_header ft_read_data ft_read_event

% the following are from the fieldtrip/external/egi directory
% TEST read_mff_header read_mff_data read_mff_event mff_getEEGFilename mff_getSummaryInfo mff_getObject mff_micros2Sample

ft_hastoolbox('egi_mff', 1);

datadir = '/home/common/matlab/fieldtrip/data/test/bug1407';
cd(datadir)

dataset = {
  'Long64ChannelWithEvents.mff'
  'NS1000sps.mff'
  'NS500Sine6Hz.mff'
  };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
headerformat = 'egi_mff_v1';
dataformat   = 'egi_mff_v1';
eventformat  = 'egi_mff_v1';

for i=1:length(dataset)
  hdr   = ft_read_header(dataset{i}, 'headerformat', headerformat);
  dat   = ft_read_data(dataset{i},   'headerformat', headerformat, 'dataformat', dataformat);
  event = ft_read_event(dataset{i},  'headerformat', headerformat, 'eventformat', eventformat);
  assert(size(dat,1)==length(hdr.label));
  assert(size(dat,2)==hdr.nSamples*hdr.nTrials);
  % keep all results to compare v1 and v2
  v1_hdr{i} = hdr;
  v1_dat{i} = dat;
  v1_evt{i} = event;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
headerformat = 'egi_mff_v2';
dataformat   = 'egi_mff_v2';
eventformat  = 'egi_mff_v2';

for i=1:length(dataset)
  hdr   = ft_read_header(dataset{i}, 'headerformat', headerformat);
  dat   = ft_read_data(dataset{i},   'headerformat', headerformat, 'dataformat', dataformat);
  event = ft_read_event(dataset{i},  'headerformat', headerformat, 'eventformat', eventformat);
  assert(size(dat,1)==length(hdr.label));
  assert(size(dat,2)==hdr.nSamples*hdr.nTrials);
  % keep all results to compare v1 and v2
  v2_hdr{i} = hdr;
  v2_dat{i} = dat;
  v2_evt{i} = event;
end

% compare the v1 results to the v2 results
for i=1:length(dataset)
  fprintf('comparing v1 and v2 for dataset %d\n', i);
  assert(isequal(v1_hdr{i}.Fs,          v2_hdr{i}.Fs)           , 'difference in hdr.Fs');
  assert(isequal(v1_hdr{i}.nChans,      v2_hdr{i}.nChans)       , 'difference in hdr.nChans');
  assert(isequal(v1_hdr{i}.nSamples,    v2_hdr{i}.nSamples)     , 'difference in hdr.nSamples');
  assert(isequal(v1_hdr{i}.nSamplesPre, v2_hdr{i}.nSamplesPre)  , 'difference in hdr.nSamplesPre');
  assert(isequal(v1_hdr{i}.nTrials,     v2_hdr{i}.nTrials)      , 'difference in hdr.nTrials');
  assert(isequal(v1_hdr{i}.label,       v2_hdr{i}.label)        , 'difference in hdr.label'); % this fails in revision 5979
  assert(isequal(v1_dat{i},  v2_dat{i})      , 'difference in data');
  assert(isequal(v1_evt{i},  v2_evt{i})      , 'difference in events'); % this fails in revision 5979
end
