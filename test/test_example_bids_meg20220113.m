function test_example_bids_meg

% MEM 4gb
% WALLTIME 00:10:00

%
%% Converting an example MEG dataset for sharing in BIDS
%
%
% The general description of how MEG data can be shared in the [BIDS format](https://bids.neuroimaging.io) is described in the paper [MEG-BIDS, the brain imaging data structure extended to magnetoencephalography](https://doi.org/10.1038/sdata.2018.110) by Guiomar Niso, et al. Scientific Data volume 5, Article number: 180110 (2018). The BIDS standard itself with all details can be found on <https://bids-standard.org>.
%
%% # CTF
%
% The CTF MEG file format consists of recordings that are represented in a directory with the extension `*.ds`. The directory contains multiple files, where the `*.res4` contains most of the header information, and the `*.meg4` files (there can be multiple) contain the data. However, also the other files contain relevant data and/or metadata, such as markers and information on how the data has been processed.
%
% In general you **should not** simply rename the directory, or mess with the files in the directory by hand, since there are multiple files that need to have a consistent file naming structure. The CTF Linux command-line software includes tools for anonimizing and renaming datasets in a consistent fashion and [this example page](/faq/how_can_i_anonymize_a_ctf_dataset/) also presents a way of doing it directly in MATLAB.
%
%% ## Example
%
% This dataset corresponds to the paper by Erik te Woerd et al. [Impaired auditory-to-motor entrainment in Parkinson’s disease](https://doi.org/10.1152/jn.00547.2016). J Neurophysiol. 2017 May 1;117(5):1853-1864. doi: 10.1152/jn.00547.2016.
%
% The complete dataset consists of MEG data recorded from 15 Parkinson patients and 15 control participants. In three of the sessions the recording was interrupted, causing those sessions to be split over 2 runs. Consequently, there are 33 CTF `*.ds` dataset directories to start with. The data from two participants was excluded from further analysis.
%
% Control subjects are coded from S101 to S115, patients are coded from S201 to S215. Besides the pseudonimized subject identifier, there is also a run number in the original file name (which is relevant for three of the recordings), the project number, the experimenter, and the date at which the recording was done. Since the date might be recognizeable for the participants and/or their family members, we will treat that as potentially identifying information and scrub it from the dataset.
%
% The data is anonimized using the `go_anonymiseDs.m` function from [this example page](/faq/how_can_i_anonymize_a_ctf_dataset/#using-matlab). This reads in the data with a low-level CTF function, scrubs potentially identifying  information, and writes it with a low-level CTF function. Some of the auxiliary files that can optionally be present (for example when the data was opened with the CTF DataEditor application) and the headlocalizer data will be ignored upon reading, and hence also not appear in the copy that is written.
%
% The anonimized version of the CTF data has a short descriptive name. Using the **[data2bids](https://github.com/fieldtrip/fieldtrip/blob/release/data2bids.m)** we copy and rename the CTF data so that it has the right directory structure and file names, and we add the metadata.
%
%
subj = {
  'S101'
  'S102'
  'S103'
  'S104'
  'S105'
  'S106'
  'S107'
  'S108'
  'S109'
  'S110'
  'S111'
  'S112'
  'S113'
  'S114'
  'S115'
  'S201'
  'S202'
  'S203'
  'S204'
  'S205'
  'S206'
  'S207'
  'S208'
  'S209'
  'S210'
  'S211'
  'S212'
  'S213'
  'S214'
  'S215'
  };


for subindx=1:numel(subj)
  for runindx=1:2
    
    % the number 3018009.04 in the dataset name refers to the DCCN project identifier
    % all data acquisition was done using study-specific participant identifiers rather than names
    % the CTF dataset name contains the recording date, but here we use a wildcard instead

    d = dir(sprintf('../original/301800904eritwoe%s00%d_1200hz_*_01.ds', subj{subindx}, runindx));
    if isempty(d)
      % for most subjects the data was recorded in a single run
      % in that case run 2 does not exist
      continue
    else
      origname = fullfile(d.folder, d.name);
      anonname = fullfile(d.folder, sprintf('%s_%d.ds', subj{subindx}, runindx));
      disp(anonname); % this is just an intermediate name, the final name will be assigned by data2bids
    end
    
    if ~isfolder(anonname)
      go_anonymiseDs(origname, anonname);
    end
    
    cfg = [];
    
    cfg.bidsroot = '../bids';  % write to the present working directory
    cfg.sub = subj{subindx};
    cfg.run = runindx;
    cfg.dataset = anonname; % this is the intermediate name
    
    cfg.datatype = 'meg';
    cfg.method = 'copy'; % the original data is in a BIDS-compliant format and can simply be copied
    
    cfg.InstitutionName             = 'Radboud University';
    cfg.InstitutionalDepartmentName = 'Donders Institute for Brain, Cognition and Behaviour';
    cfg.InstitutionAddress          = 'Kapittelweg 29, 6525 EN, Nijmegen, The Netherlands';
    
    % required for dataset_description.json
    cfg.dataset_description.Name                = 'Impaired auditory-to-motor entrainment in Parkinson''s disease';
    cfg.dataset_description.BIDSVersion         = 'v1.5.0';
    
    % optional for dataset_description.json
    cfg.dataset_description.License             = 'RU-DI-HD-1.0';
    cfg.dataset_description.Authors             = 'Te Woerd ES, Oostenveld R, de Lange FP, Praamstra P';
    cfg.dataset_description.ReferencesAndLinks  = {'https://doi.org/10.1152/jn.00547.2016'};
    cfg.dataset_description.EthicsApprovals     = 'DCCN 3018009.04';
    
    cfg.TaskDescription = 'The experiment consisted of an auditory target detection task. Three cue tones of different frequencies (900, 1,100 and 1,300 Hz; 50-ms duration) predicted the probability (10, 30, 50%) of the next stimulus being the target tone (2,000 Hz; 50-ms duration) (Fig. 1A). The pitch of standard tones was chosen to be separated by 200 Hz and the target differed 700 Hz from the highest standard tone, which enabled subjects to detect the target tone easily. Participants were informed about the meaning of the standard tones and instructed to press the response button as swift as possible with the index finger of their dominant hand.';
    cfg.task = 'auditory';
    
    try
      data2bids(cfg);
    catch
      % this is probably because the output dataset already exists
      % this is due to running the script multiple times
      disp(lasterr)
    end
    
    
  end % for each run
end % for each subject
%
% Following this conversion, the `participant.tsv` file only contains the subject identifiers. There is more information in an Excel file, such as the age, sex, clinical state, questionaire data, lab notes. In principle it would be possible to copy that into the MATLAB script and use data2bids to export that, but in this case it was easier to use "save-as" and then select "tab-delimited-text" as the output format.
