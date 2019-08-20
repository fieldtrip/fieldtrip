function test_issue1196

% MEM 2gb
% WALLTIME 00:20:00
% DEPENDENCY data2bids

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1196'));

if isfolder('/home/common/matlab/fieldtrip/data/test/issue1196/bids')
  system('rm -rf /home/common/matlab/fieldtrip/data/test/issue1196/bids');
end

%%

emgfile = 'raw/006_3013065.02_rest1.vhdr';
t1wfile = 'raw/sub-139075/015-T1w_MPR/00175_1.3.12.2.1107.5.2.43.66068.2019080415440698293169043.IMA'; % any slice is fine
t2wfile = 'raw/sub-139075/017-T2w_SPC/00011_1.3.12.2.1107.5.2.43.66068.2019080415500178296969403.IMA';

%%
% the first part is general to the different recordings

cfg = [];
cfg.bidsroot = 'bids';
cfg.sub = '139075';
cfg.ses = 'mri'; % all data was recorded in a single MRI session

cfg.InstitutionName             = 'Radboud University';
cfg.InstitutionalDepartmentName = 'Donders Institute for Brain, Cognition and Behaviour';
cfg.InstitutionAddress          = 'Kapittelweg 29, 6525 EN, Nijmegen, The Netherlands';

% REQUIRED
cfg.dataset_description.Name                = 'n/a';
cfg.dataset_description.BIDSVersion         = '1.2';
% OPTIONAL
cfg.dataset_description.License             = 'n/a';
cfg.dataset_description.Authors             = 'n/a';
cfg.dataset_description.Acknowledgements    = 'n/a';
cfg.dataset_description.Funding             = 'n/a';
cfg.dataset_description.ReferencesAndLinks  = 'n/a';
cfg.dataset_description.DatasetDOI          = 'n/a';

%%
% convert the T1w dicom files

cfg.dataset = t1wfile;
cfg.datatype = 'T1w';
data2bids(cfg);

%%
% convert the T2w dicom files

cfg.dataset = t2wfile;
cfg.datatype = 'T2w';
data2bids(cfg);

%%
% copy the EMG

% these are general fields and need to be added for EMG
% for MRI they can (hopefully) be read from the DICOM headers
cfg.Manufacturers = 'BrainProducts';
cfg.ManufacturersModelName = 'BrainAmp MR plus';
cfg.DeviceSerialNumber = 'n/a';
cfg.SoftwareVersions = 'n/a';

% FIXME we should agree upon the metadata fields for the JSON
% cfg.emg.SamplingFrequency         % this will be determined from the file
% cfg.emg.RecordingDuration         % this will be determined from the file, in seconds
% cfg.emg.RecordingType             % this will probably be determined from the file, continuous
cfg.emg.PowerLineFrequency = 50;
% cfg.emg.HardwareFilters
% cfg.emg.SoftwareFilters
cfg.emg.EMGChannelCount = 2;
cfg.emg.EOGChannelCount = pi;
cfg.emg.ECGChannelCount = 42;
cfg.emg.ElectrodeManufacturer = 'Kendall';
cfg.emg.ElectrodeManufacturersModelName = 'schuim, 24mm';
% cfg.emg.EMGPlacementScheme
% cfg.emg.EMGReference
% cfg.emg.EMGGround

cfg.method = 'copy'; % no need to convert
cfg.dataset = emgfile;
cfg.datatype = 'emg';

data2bids(cfg);


