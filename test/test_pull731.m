function test_pull731

% MEM 12gb
% WALLTIME 01:00:00
% DEPENDENCY ft_read_header ft_read_data ft_read_event

datadir = dccnpath('/home/common/matlab/fieldtrip/data/test/pull731');

dataset = {
  fullfile(datadir, './OtherFilesRobert/bug1427/Long64ChannelWithEvents.mff')
  fullfile(datadir, './OtherFilesRobert/bug1427/NS1000sps.mff')
  fullfile(datadir, './OtherFilesRobert/bug1427/NS500Sine6Hz.mff')
  fullfile(datadir, './OtherFilesRobert/bug629/pilot05_test 20110120 1433.mff')
  fullfile(datadir, './Bugs/SPBI023 20150414 1357.mff')
  fullfile(datadir, './Continuous with Video/ZZ2_018 0416 2110.mff')
  fullfile(datadir, './Grand Average Multiple Subjects/GNG_F_Day_1_GAve.mff')
  fullfile(datadir, './Individual Averaging multiple categories/LLL_01.1_T108_0691.ave.mff')
  fullfile(datadir, './Individual Averaging multiple subjects/GNG_F_Day_1_Combined_blc_ref.mff')
  fullfile(datadir, './OtherFilesArnaud/9999_20160309_011903.mff')
  fullfile(datadir, './OtherFilesArnaud/Nicolas_CerCo_SI_STE_test1_20151009_034300.mff')
  fullfile(datadir, './OtherFilesArnaud/coma_LA_MOH1_20160411_060749.mff')
  fullfile(datadir, './OtherFilesArnaud/data.mff')
  fullfile(datadir, './OtherFilesArnaud/multiSubj_seg_LLL_06.cmbave.mff')
  fullfile(datadir, './OtherFilesDavid/01_024 0531 1145_seg_fil_bcr_ave_WITH_AUTONOMOUS.mff')
  fullfile(datadir, './OtherFilesDavid/4ms_5uV.nsr.mff')
  fullfile(datadir, './OtherFilesDavid/ATD256_1.mff')
  fullfile(datadir, './OtherFilesDavid/GNG2_014_1s_cln.seg.mff')
  fullfile(datadir, './OtherFilesDavid/VTD_7Ss_bcr.gav.mff')
  fullfile(datadir, './OtherFilesDavid/VTD_993.1.ses .mff')
  fullfile(datadir, './OtherFilesRobert/original/eeg/egi/NS500Sine6Hz.mff')
  fullfile(datadir, './Segmented with multiple categories/NIA_P_013.ses.fil.seg.mff')
  fullfile(datadir, './Treys_files/MMI_HC1_20180314_093330_physio_only.mff')
  fullfile(datadir, './Treys_files/MMVTD_Continuous_EEG.mff')
  fullfile(datadir, './Unprocessed Continuous/128 channels/GNG2_002_1v_cln.nsf.mff')
  fullfile(datadir, './Unprocessed Continuous/256 channels/VTD_993.1.ses.mff')
  fullfile(datadir, './Unprocessed Continuous/32 channels/NIA_333ms_HCGSN32_test01.mff')
  fullfile(datadir, './Unprocessed Continuous/32 channels/test.mff')
  fullfile(datadir, './Unprocessed Continuous/64 channels/NIA_333msHCGSN64_test01.mff')
  };


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
headerformat = 'egi_mff_v1';
dataformat   = 'egi_mff_v1';
eventformat  = 'egi_mff_v1';

for i=1:length(dataset)
  try
    hdr   = ft_read_header(dataset{i}, 'headerformat', headerformat);
    dat   = ft_read_data(dataset{i},   'headerformat', headerformat, 'dataformat', dataformat);
    event = ft_read_event(dataset{i},  'headerformat', headerformat, 'eventformat', eventformat);
    assert(size(dat,1)==length(hdr.label));
    assert(size(dat,2)==hdr.nSamples*hdr.nTrials);
    % keep all results to compare v1 and v2
    v1_hdr{i} = hdr;
    v1_dat{i} = dat;
    v1_evt{i} = event;
  catch
    warning('the egi_mff_v1 implementation fails for %s', dataset{i});
    % there are quite a few that the egi_mff_v2 implementation fails to read
    v1_hdr{i} = [];
    v1_dat{i} = [];
    v1_evt{i} = [];
  end
end

% the purpose for the following code is try ensure that the v1 implementation does not break
% these numbers are hardcoded for the datasets defined above
assert(v1_hdr{1}.nSamples*v1_hdr{1}.nTrials==3087,  'incorrect number of samples')
assert(v1_hdr{2}.nSamples*v1_hdr{2}.nTrials==60660, 'incorrect number of samples')
assert(v1_hdr{3}.nSamples*v1_hdr{3}.nTrials==77113, 'incorrect number of samples')
assert(v1_hdr{1}.nChans==65,   'incorrect number of samples')
assert(v1_hdr{2}.nChans==257,  'incorrect number of samples')
assert(v1_hdr{3}.nChans==257,  'incorrect number of samples')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
headerformat = 'egi_mff_v2';
dataformat   = 'egi_mff_v2';
eventformat  = 'egi_mff_v2';

for i=1:length(dataset)
  try
    hdr   = ft_read_header(dataset{i}, 'headerformat', headerformat);
    dat   = ft_read_data(dataset{i},   'headerformat', headerformat, 'dataformat', dataformat);
    event = ft_read_event(dataset{i},  'headerformat', headerformat, 'eventformat', eventformat);
    assert(size(dat,1)==length(hdr.label));
    assert(size(dat,2)==hdr.nSamples*hdr.nTrials);
    % keep all results to compare v1 and v2
    v2_hdr{i} = hdr;
    v2_dat{i} = dat;
    v2_evt{i} = event;
  catch
    warning('the egi_mff_v2 implementation fails for %s', dataset{i});
    % there are quite a few that the egi_mff_v2 implementation fails to read
    v2_hdr{i} = [];
    v2_dat{i} = [];
    v2_evt{i} = [];
  end
end

%% compare the v1 results to the v2 results (where possible)
for i=1:length(dataset)
  if isempty(v1_hdr{i}) || isempty(v2_hdr{i})
    fprintf('not comparing v1 and v2 for dataset %d\n', i);
  else
    fprintf('comparing v1 and v2 for dataset %d\n', i);
    assert(isequal(v1_hdr{i}.Fs,          v2_hdr{i}.Fs)           , 'difference in hdr.Fs');
    assert(isequal(v1_hdr{i}.nChans,      v2_hdr{i}.nChans)       , 'difference in hdr.nChans');
    assert(isequal(v1_hdr{i}.nSamples,    v2_hdr{i}.nSamples)     , 'difference in hdr.nSamples');
    assert(isequal(v1_hdr{i}.nSamplesPre, v2_hdr{i}.nSamplesPre)  , 'difference in hdr.nSamplesPre');
    assert(isequal(v1_hdr{i}.nTrials,     v2_hdr{i}.nTrials)      , 'difference in hdr.nTrials');
    assert(isequal(v1_hdr{i}.label,       v2_hdr{i}.label)        , 'difference in hdr.label'); % this fails in revision 5979
    assert(isequal(v1_dat{i},             v2_dat{i})              , 'difference in data');
    assert(length(v1_evt{i})==length(v2_evt{i}) , 'difference in number of events');
    assert(isequal(sort([v1_evt{i}.sample]), sort([v2_evt{i}.sample])), 'difference in event.sample');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
headerformat = 'egi_mff_v3';
dataformat   = 'egi_mff_v3';
eventformat  = 'egi_mff_v3';

for i=1:length(dataset)
  try
    hdr   = ft_read_header(dataset{i}, 'headerformat', headerformat);
    dat   = ft_read_data(dataset{i},   'headerformat', headerformat, 'dataformat', dataformat);
    event = ft_read_event(dataset{i},  'headerformat', headerformat, 'eventformat', eventformat);
    assert(size(dat,1)==length(hdr.label));
    assert(size(dat,2)==hdr.nSamples*hdr.nTrials);
    % keep all results to compare v1 and v3
    v3_hdr{i} = hdr;
    v3_dat{i} = dat;
    v3_evt{i} = event;
  catch
    warning('the egi_mff_v3 implementation fails for %s', dataset{i});
    v3_hdr{i} = [];
    v3_dat{i} = [];
    v3_evt{i} = [];
  end
end

%% compare the v1 results to the v3 results (where possible)
for i=1:length(dataset)
  if isempty(v1_hdr{i}) || isempty(v3_hdr{i})
    fprintf('not comparing v1 and v3 for dataset %d\n', i);
  else
    % there are quite some differences
    % there is no point in giving errors, so use wassert instead
    fprintf('comparing v1 and v3 for dataset %d\n', i);
    wassert(isequal(v1_hdr{i}.Fs,          v3_hdr{i}.Fs)           , 'difference in hdr.Fs');
    wassert(isequal(v1_hdr{i}.nChans,      v3_hdr{i}.nChans)       , 'difference in hdr.nChans');
    wassert(isequal(v1_hdr{i}.nSamples,    v3_hdr{i}.nSamples)     , 'difference in hdr.nSamples');
    wassert(isequal(v1_hdr{i}.nSamplesPre, v3_hdr{i}.nSamplesPre)  , 'difference in hdr.nSamplesPre');
    wassert(isequal(v1_hdr{i}.nTrials,     v3_hdr{i}.nTrials)      , 'difference in hdr.nTrials');
    wassert(isequal(v1_hdr{i}.label,       v3_hdr{i}.label)        , 'difference in hdr.label'); % this fails in revision 5979
    wassert(isequal(v1_dat{i},             v3_dat{i})              , 'difference in data');
    wassert(length(v1_evt{i})==length(v3_evt{i}) , 'difference in number of events');
    wassert(isequal(sort([v1_evt{i}.sample]), sort([v3_evt{i}.sample])), 'difference in event.sample');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wassert(condition, message)
if ~condition
  warning(message)
end
