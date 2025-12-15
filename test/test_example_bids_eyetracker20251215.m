function test_example_bids_eyetracker_20251215

% MEM 2gb
% WALLTIME 00:20:00
% DATA public

% The [BIDS standard](https://www.bids-standard.org) does currently not specify how to represent eye tracker data, but an extension is discussed [here](https://github.com/bids-standard/bids-specification/discussions/2218). This example - and the support that is implemented in the **[data2bids](/reference/data2bids)** function - should be considered as a preliminary proposal to help researchers with their existing data.
%
% Eye tracking data can be stored in the BIDS representation similar to [physiological data](https://bids-specification.readthedocs.io/en/latest/modality-specific-files/physiological-recordings.html). This not only includes continuously sampled variables like the gaze position and pupil diameter to be stored in the physio.tsv file, but also allows for presentation (stimulus and response events) and saccades to be represented in the physioevents.tsv file.
%
% The physio.json file should contain additional metadata information such as whether both eyes or only one was tracked, the sapling rate, if and how pupil diameter is quantified, what calibration process was used, how to interpret the gaze position (pixels, degrees), whether the origin is at the center of the screen or the upper left corner, etc.
%
% All data for the following examples is available from our [download server](https://download.fieldtriptoolbox.org/example/bids_eyetracker/).
%
% The examples includes the original log files under the “original” directory. It also has a copy of the script to do the conversion under “code”. The reorganized data is under the “bids” directory. According to [the documentation](https://bids-specification.readthedocs.io/en/stable/02-common-principles.html#source-vs-raw-vs-derived-data) the original data can be added to the BIDS dataset under the “sourcedata” directory and code can be added to the “code” directory. This way no information is lost and the conversion/reorganization is fully reproducible.
%
% In general: if your original data is in a different format than the BIDS representation (e.g., DICOM instead of NIfTI, or the original binary data from the eye tracker software), you probably want  to keep a copy of the original data as well, e.g., on a data acquisition collection on the Donders Repository. If it is in the same format (which is not the case here), you could simply delete the original data after conversion. In either case - your own analyses and the shared data would be based on the BIDS representation.
%
% There are numerous eye tracking systems, such as SR Research EyeLink, Tobii, EyeTech, GazePoint, SmartEye and Pupil Labs. We cannot provide examples for all of them, but in principle data from each system would be stored in the same standard BIDS representation, allowing it to be reused by others.
%
% There is a draft BIDS extension proposal (BEP) to extend the specification with eye tracker data, you can find the discussion on [google docs](https://docs.google.com/document/d/1eggzTCzSHG3AEKhtnEDbcdk-2avXN6I94X8aUPEBVsw/edit#).
%
%% # SR Research EyeLink
%
% The SR Research  EyeLink system records the data in a binary file with the extension `.edf`. This is not to be confused with the European Data Format (https://www.edfplus.info) for EEG and other biological and physical signals. The proprietary EyeLink edf format cannot be read in open source software, but SR Research provides a tool called `edf2asc.exe` that converts the data to ASCII format. The ASCII format is still difficult to work with and represents a mixture of continuous gaze and pupil diameter parameters with discontinuous events, such as fixation, loss of tracking, etc.
%
% FieldTrip can read the EyeLink `.asc` format. The gaze and/or pupil diameter are represented as continuous raw channels (accessible using **[ft_read_data](/reference/fileio/ft_read_data)** and **[ft_preprocessing](/reference/ft_preprocessing)**), the other information is represented as events (using **[ft_read_event](/reference/fileio/ft_read_event)** and **[ft_definetrial](/reference/ft_definetrial)**). For more details you can see [this page](/getting_started/eyetracker/eyelink). In the examples below we use **[data2bids](/reference/data2bids)** to read the `.asc` files and convert the data to a simple `_eyetracker.tsv` file for the continuous gaze and pupil diameter, and an `_events.tsv` file for the discontinuous events.
%
%% ## Short example
%
% The data for this example is available from our [download server](https://download.fieldtriptoolbox.org/example/bids_eyetracker/eyelink_short).
%
%% this is an example with a short calibration recording

% the same data is converted multiple times to demonstrate the directory layout and the participants

cfg = [];

cfg.InstitutionName             = 'Radboud University';
cfg.InstitutionalDepartmentName = 'Donders Institute for Brain, Cognition and Behaviour';
cfg.InstitutionAddress          = 'Kapittelweg 29, 6525 EN, Nijmegen, The Netherlands';

% required for dataset_description.json
cfg.dataset_description.Name        = 'Eyelink example';
cfg.dataset_description.BIDSVersion = 'unofficial extension';

% optional for dataset_description.json
cfg.dataset_description.License             = 'n/a';
cfg.dataset_description.Authors             = 'n/a';
cfg.dataset_description.Acknowledgements    = 'n/a';
cfg.dataset_description.Funding             = 'n/a';
cfg.dataset_description.ReferencesAndLinks  = 'n/a';
cfg.dataset_description.DatasetDOI          = 'n/a';

cfg.dataset   = dccnpath(fullfile('/project/3031000.02/external/download/example/bids_eyetracker/eyelink_short/original', 'ashcal.asc'));
cfg.method    = 'convert'; % the eyelink EDF format is not supported, convert the ASC to plain TSV
cfg.bidsroot  = fullfile(tempdir, 'bids');
cfg.suffix    = 'physio';
cfg.task      = 'calibration';

% this is general metadata that ends up in the _physio.json file
cfg.PhysioType            = 'eyetrack';
cfg.TaskDescription       = 'Short calibration procedure';
cfg.Manufacturer          = 'SR Research';
cfg.ManufacturerModelName = 'Eyelink 1000';

% convert the data from the first (and only one) subject
cfg.sub = '01';
cfg.participants.age = 23;
cfg.participants.sex = 'M';
data2bids(cfg);

% convert the same data again, pretending that it is subject 2
cfg.sub = '02';
cfg.participants.age = 25;
cfg.participants.sex = 'F';
data2bids(cfg);

% convert the same data again, pretending that it is subject 3
cfg.sub = '03';
cfg.participants.age = 31;
cfg.participants.sex = 'M';
data2bids(cfg);
%
%% ## Long example
%
% The data for this example is available from our [download server](https://download.fieldtriptoolbox.org/example/bids_eyetracker/eyelink_long).
%
%% this section specifies the datafiles

% The raw sourcedata for this example can be downloaded from
% https://doi.org/10.34973/3mdw-jn48 as Experiment 2.

filename = {
  '1_TO_1.asc'
  '1_TO_2.asc'
  '2_FH_1.asc'
  '2_fh_2.asc'
  '3_ac_1.asc'
  '3_ac_2.asc'
  % '4_ym_1.asc'
  % '4_ym_2.asc'
  % '5_kj_1.asc'
  % '5_kj_2.asc'
  % '6_BS_1.asc'
  % '6_BS_2.asc'
  % '7_aa_1.asc'
  % '7_aa_2.asc'
  % '8_SG_1.asc'
  % '8_sg_2.asc'
  % '9_TE_1.asc'
  % '9_TE_2.asc'
  % '9_te_11.asc'
  % '10_es_1.asc'
  % '10_es_2.asc'
  % '11_AR_1.asc'
  % '11_AR_2.asc'
  % '12_LK_1.asc'
  % '12_LK_11.asc'
  % '12_LK_2.asc'
  % '12_LK_22.asc'
  % '13_JL_1.asc'
  % '13_JL_2.asc'
  % '14_vs_1.asc'
  % '14_vs_2.asc'
  % '15_CL_1.asc'
  % '15_CL_2.asc'
  % '16_MN_1.asc'
  % '16_MN_2.asc'
  % '17_NT_1.asc'
  % '17_NT_2.asc'
  % '18_CL_11.asc'
  % '18_CL_2.asc'
  % '19_AE_1.asc'
  % '19_AE_2.asc'
  % '20_AK_2.asc'
  % '20_PK_1.asc'
  % '21_SY_1.asc'
  % '21_SY_2.asc'
  % '22_VG_1.asc'
  % '22_VG_21.asc'
  % '22_VG_22.asc'
  % '23_DK_2.asc'
  % '23_dk_1.asc'
  % '24_TS_1.asc'
  % '24_TS_2.asc'
  };

% The filenames have a sequence number and a two-letter subject identifier in them. The case of the
% identifier is not totally consistent. Furthermore, some subjects have two, some
% have three recordings.

sourcepath = dccnpath('/project/3031000.02/external/download/example/bids_eyetracker/eyelink_long/original/ascData');
targetpath = fullfile(tempdir, 'bids');

%% find the unique subject identifiers

subjid = {};
for i=1:numel(filename)
  part = split(filename{i}, '_');
  subjid{end+1} = upper(part{2});
end
subjid = unique(subjid, 'stable'); % keep them in the original order

%% loop over subjects and convert the datafiles for each subject

cfg = [];

cfg.InstitutionName             = 'Radboud University';
cfg.InstitutionalDepartmentName = 'Donders Institute for Brain, Cognition and Behaviour';
cfg.InstitutionAddress          = 'Kapittelweg 29, 6525 EN, Nijmegen, The Netherlands';

% required for dataset_description.json
cfg.dataset_description.Name                = 'Visual Stability: predictive remapping of features beyond saccadic targets';
cfg.dataset_description.BIDSVersion         = 'unofficial extension';

% optional for dataset_description.json
cfg.dataset_description.Authors             = 'Tao He, Matthias Fritsche, Floris de Lange';
cfg.dataset_description.DatasetDOI          = 'http://hdl.handle.net/11633/di.dccn.DSC_3018034.01_694';
cfg.dataset_description.License             = 'RU-DI-HD-1.0';
cfg.dataset_description.Acknowledgements    = 'n/a';
cfg.dataset_description.Funding             = 'n/a';
cfg.dataset_description.ReferencesAndLinks  = 'n/a';

for i=1:numel(subjid)
  fileselection = find(contains(filename, ['_'  subjid{i} '_'], 'IgnoreCase', true));
  
  for j=1:numel(fileselection)
    [p, f, x] = fileparts(filename{fileselection(j)});  % split the path, filename and extension
    part = split(f, '_');                               % split the filename
    run = str2double(part{3});
    
    % this is the data in ASCII format
    cfg.dataset = fullfile(sourcepath, [f x]);
    cfg.method = 'convert';
    
    % these are used to construct the directory and file name
    cfg.bidsroot = targetpath;
    cfg.sub = subjid{i};
    cfg.run = run;
    cfg.task = 'adaptation';
    cfg.suffix = 'physio';

    % this ends up in the physio.json file
    cfg.physio.PhysioType = 'eyetrack';
    
    % this is additional information that ends up in the sidecar JSON file
    cfg.TaskDescription = 'orientation adaptation paradigm';
    
    data2bids(cfg);
    
  end % for fileselection
end % for subjid
%
%% # SMI
%
% The SMI eye tracker stores the raw data in an `.idf` file. That file cannot be read (easily) with other analysis software, hence the SMI specific *IDF converter* program from the iView tools needs to be used to convert it to ascii. Please note that this is Windows software that comes with the eye tracker. At the Donders Centre for Cognitive Neuroimaging (DCCN) in Nijmegen you can ask Paul for more details on the SMI eye trackers.
%
%% ## Example
%
% In the following example we are converting two runs of eye tracker data for two subjects. The data was recorded at the DCCN. Since additional information is missing (e.g., units, origin, calibration procedure), the metadata is very sparse.
%
% The data for this example is available from our [download server](https://download.fieldtriptoolbox.org/example/bids_eyetracker/smi).
%
datadir = dccnpath('/project/3031000.02/external/download/example/bids_eyetracker/smi/original');
filename = {
  fullfile(datadir, 'pp23671_rest1_samples.txt')
  fullfile(datadir, 'pp23671_task1_samples.txt')
  fullfile(datadir, 'pp31237_rest1_samples.txt')
  fullfile(datadir, 'pp31237_task1_samples.txt')
  };

% note that the original filename includes the subject identifier and the task
% these are used further down and included in the output

cfg = [];

cfg.InstitutionName             = 'Radboud University';
cfg.InstitutionalDepartmentName = 'Donders Institute for Brain, Cognition and Behaviour';
cfg.InstitutionAddress          = 'Kapittelweg 29, 6525 EN, Nijmegen, The Netherlands';

% required for dataset_description.json
cfg.dataset_description.Name        = 'SMI example';
cfg.dataset_description.BIDSVersion = 'unofficial extension';

% optional for dataset_description.json
cfg.dataset_description.License             = 'n/a';
cfg.dataset_description.Authors             = 'n/a';
cfg.dataset_description.Acknowledgements    = 'n/a';
cfg.dataset_description.Funding             = 'n/a';
cfg.dataset_description.ReferencesAndLinks  = 'n/a';
cfg.dataset_description.DatasetDOI          = 'n/a';

cfg.Manufacturer          = 'SMI';
cfg.ManufacturerModelName = 'iView X MRI-LR';

cfg.method    = 'convert'; % the SMI-specific format is not supported, convert it to plain TSV
cfg.bidsroot  = fullfile(tempdir, 'bids');
cfg.suffix    = 'physio';

for i=1:4
  cfg.dataset   = filename{i};
  
  % split the filename to get the subject identifier and the task
  [p, f, x] = fileparts(filename{i});
  piece = split(f, '_');
  
  cfg.sub   = piece{1};
  cfg.task  = piece{2};
  
  data2bids(cfg);
end
%
%% # TOBII
%
% We currently do not have a complete example for converting TOBII eye tracker data to BIDS, but have collected some information that is posted here for reference.
%
% The data from the TOBII eye tracker can be exported in `.tsv` or in `.xlsx` format. Moreover, the data can be exported to either have a file for each subject, or with all the subjects in one file. The TOBI studio software allows to add additional information in the exported files.
