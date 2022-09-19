function test_example_bids_emg

% MEM 4gb
% WALLTIME 00:20:00

%% Converting an example EMG dataset for sharing in BIDS
%
% The [BIDS standard](https://bids.neuroimaging.io) does currently not specify how to represent EMG recordings. This example - and the support that is implemented in the **[data2bids](https://github.com/fieldtrip/fieldtrip/blob/release/data2bids.m)** function - should be considered as a preliminary proposal to help researchers with their existing data. This example may also serve to start a discussion on whether and how this data type should be added to the [BIDS specification](http://bids-specification.readthedocs.io/).
%
% EMG is often recorded along with other experiments, such as MEG or functional MRI recordings, but can also be recorded by itself. EMG can for example be used to check compliance with the task, to detect movements and compute reaction times, or to check the activation pattern of specific muscles in a task.
%
% EMG is not part as an individual data type in the BIDS specification. In the **[data2bids](https://github.com/fieldtrip/fieldtrip/blob/release/data2bids.m)** function we are mostly treating it similar to EEG data, with only some differences in the metadata. E.g., the |EEGPlacementScheme|, |CapManufacturer|, and |CapManufacturersModelName| do not apply.  Other than that, we have implemented it the same, with the same data formats and the same way to code events.
%
%% # Example
%
% The example data is available from our [FTP server](ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/example/bids_emg/)
%
% The example includes the original log files under the “original” directory. It also has a copy of the script to do the conversion under “code”. The reorganized data is under the “bids” directory. According to [the documentation](https://bids-specification.readthedocs.io/en/stable/02-common-principles.html#source-vs-raw-vs-derived-data) the original data can be added to the BIDS dataset under the “sourcedata” directory and code can be added to the “code” directory. This way no information is lost and the conversion/reorganization is fully reproducible.
%
% In general: if your original data is in a different format than the BIDS representation (e.g., DICOMs instead of NIfTI, or data from an EMG system in an unsupported format), you might want to keep the sourcedata, e.g., on a data acquisition collection on the Donders Repository. If it is in the same format (as here), you could simply delete the copy in sourcedata. Your own analyses and the sharing of the data would be based on the BIDS representation.
%
cfg = [];

cfg.InstitutionName             = 'Radboud University';
cfg.InstitutionalDepartmentName = 'Donders Institute for Brain, Cognition and Behaviour';
cfg.InstitutionAddress          = 'Kapittelweg 29, 6525 EN, Nijmegen, The Netherlands';

% required for dataset_description.json
cfg.dataset_description.Name                = 'EMG example';
cfg.dataset_description.BIDSVersion         = 'unofficial extension';

% optional for dataset_description.json
cfg.dataset_description.License             = 'n/a';
cfg.dataset_description.Authors             = 'n/a';
cfg.dataset_description.Acknowledgements    = 'n/a';
cfg.dataset_description.Funding             = 'n/a';
cfg.dataset_description.ReferencesAndLinks  = 'n/a';
cfg.dataset_description.DatasetDOI          = 'n/a';

cfg.method = 'copy'; % the original data is in a BIDS -compliant format and can simply be copied

cfg.bidsroot = fullfile(tempdir, 'bids');  % write to the some scratch directory
cfg.datatype = 'emg';
cfg.sub = 'NT05';

% these are general fields
cfg.Manufacturer           = 'BrainProducts';
cfg.ManufacturersModelName = 'BrainAmp MR plus';

% these are EMG specific fields
cfg.emg.EMGPlacementScheme = 'electrode pairs were placed over various muscles';
cfg.emg.EMGReference       = 'bipolar';
cfg.emg.EMGGround          = 'unclear where the ground electrode was placed on the body';
cfg.emg.EMGChannelCount    = 13; % it might actually be that some channels are recording something else

% there are two recordings of the EMG, during the RS scan and during the GABA scan
% in both EMG recordings the task for the subject was the same

cfg.TaskDescription = 'The subject was lying in the MRI scanner during a resting-state BOLD scan, while EMG was being recorded';
cfg.task = 'resting';
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/example/bids_emg/original/nt05_EMG_SESS1_RS_RAW.vhdr');
cfg.run = 1;

% write the first EMG recording as run-1
data2bids(cfg);

cfg.TaskDescription = 'The subject was lying in the MRI scanner during a GABA MRS scan, while EMG was being recorded';
cfg.task = 'resting';
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/example/bids_emg/original/nt05_EMG_SESS1_GABA_RAW.vhdr');
cfg.run = 2;

% write the first EMG recording as run-2
data2bids(cfg);
