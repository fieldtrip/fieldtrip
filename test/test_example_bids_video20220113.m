function test_example_bids_video

% MEM 4gb
% WALLTIME 00:10:00

%
%% Converting an example video dataset for sharing in BIDS
%
% The [BIDS standard](https://bids.neuroimaging.io) does currently not specify how to represent video recordings. This example - and the support that is implemented in the **[data2bids](https://github.com/fieldtrip/fieldtrip/blob/release/data2bids.m)** function - should be considered as a preliminary proposal to help researchers with their existing data. This example may also serve to start a discussion on whether and how this data type should be added to the [BIDS specification](http://bids-specification.readthedocs.io/).
%
% In many experiments video data can be recorded along with other measurements. This is often used to off-line screen the subjects behavior and to annotate the video, i.e. note the time segments with specific (desired or undesired) behavior. Neither the video recording  (which can be considered as raw data), nor the annotations (which can be considered as derived data) are currently formally part of BIDS, although it is not so difficult to come up with a way that these can be represented in a BIDS-like fashion.
%
%% # Example
%
% The example data is available from our [FTP server](ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/example/bids_video/).
%
% The example includes the original video recording under the “original” directory. It also has a copy of the script to do the conversion under “code”. The reorganized data is under the “bids” directory. According to [the documentation](https://bids-specification.readthedocs.io/en/stable/02-common-principles.html#source-vs-raw-vs-derived-data) the original data can be added to the BIDS dataset under the “sourcedata” directory and code can be added to the “code” directory. This way no information is lost and the conversion/reorganization is fully reproducible.
%
% In general: if your original data is in a different format than the BIDS representation (e.g., DICOM instead of NIFTI), you probably want to keep a copy of the original data, e.g., on a data acquisition collection on the Donders Repository. If it is in the same format like here - since we are not converting the video file, you could simply delete the original files after conversion. In either case - your own analyses and the shared data would be based on the BIDS representation.
%
cfg = [];

cfg.InstitutionName             = 'Radboud University';
cfg.InstitutionalDepartmentName = 'Donders Institute for Brain, Cognition and Behaviour';
cfg.InstitutionAddress          = 'Kapittelweg 29, 6525 EN, Nijmegen, The Netherlands';

% this is required for dataset_description.json
cfg.dataset_description.Name                = 'Video example';
cfg.dataset_description.BIDSVersion         = 'unofficial extension';

% this is optional for dataset_description.json
cfg.dataset_description.License             = 'n/a';
cfg.dataset_description.Authors             = 'n/a';
cfg.dataset_description.Acknowledgements    = 'n/a';
cfg.dataset_description.Funding             = 'n/a';
cfg.dataset_description.ReferencesAndLinks  = 'n/a';
cfg.dataset_description.DatasetDOI          = 'n/a';

% provide some metadata about the task that was performed
cfg.TaskDescription = 'The subject was walking, stopping, turning, etcetera, while wearing a mobile EEG system on his head and a mobile NIRS system in his pocket';

cfg.method    = 'copy'; % video should simply be copied, not converted
cfg.dataset   = dccnpath('/home/common/matlab/fieldtrip/data/ftp/example/bids_video/original/Video_20190411_210224.mp4');
cfg.bidsroot  = fullfile(tempdir, 'bids');
cfg.datatype  = 'video';
cfg.sub       = '01';
cfg.task      = 'walking';

% construct an empty table that has the required columns 'onset' and 'duration' (in seconds)
% this uses the optional 'trial_type' column to indicate what happened
onset      = nan;
duration   = nan;
trial_type = {''};
cfg.events = table(onset, duration, trial_type);

% add the actual  events with onset, duration and trial_type
cfg.events(1,:) = {12, 2, 'standing up'};
cfg.events(2,:) = {15, 11, 'walking'};
cfg.events(3,:) = {26, 3, 'right turn'};
cfg.events(4,:) = {29, 10, 'walking'};

data2bids(cfg);
%
% This example also shows how to add annotations of the video to the `_events.tsv` file. If you do not have annotations (yet), you would leave the `cfg.events` field empty.
%
% In case you would have multiple subsequent videos in a single session, you would use `cfg.run=1, 2, ...` to specify the successive recordings.
