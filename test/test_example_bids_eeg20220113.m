function test_example_bids_eeg

% MEM 4gb
% WALLTIME 00:10:00



%
%% Converting an example EEG dataset for sharing in BIDS
%
% This example describes how to use **[data2bids](https://github.com/fieldtrip/fieldtrip/blob/release/data2bids.m)** to convert an EEG dataset for sharing according to the [BIDS standard](https://bids.neuroimaging.io).
%
% If you want to share data, there are multiple things to consider. For example the file format of the data, the place to share the data (openneuro/zenodo/figshare/dataverse), the access mechanism for the data (ftp/http/sftp), the license or data use agreement, whether all data or only part of the data is shared, using pseudonyms in the description of the data, scrubbing the date and time of recording, removing identifying features from the data, etc.
%
% In this example we will only be dealing with the format in which the data is organized (over directories) and stored (in files), for which we use BIDS. Please see the [BIDS website](https://bids.neuroimaging.io) for background information and the [BIDS specification](https://bids-specification.readthedocs.io/en/stable/) for further details on the file and directory organization.
%
% We will describe two approaches: in the first the data files are is kept in their original format, in the second approaches the files are explicitly converted to BrainVision (the recommended format for EEG data in BIDS). Prior to conversion the data (in principle) comprises 10 files, one EEG file per subject. After conversion there are 52 or 72 files for the two approaches, which includes the EEG data and the sidecar files with metadata.
%
% It is important that you use appropriate tools. BIDS stores metadata in TSV and JSON files, which makes then not only machine readable, but also human readable. A good graphical text/code editor helps you to navigate through the full directory structure and check or update the content of individual files. We recommend the [atom editor](http://atom.io/), but there are good [alternatives](https://alternativeto.net/software/atom/).
%
% The BIDS standard for EEG is introduced in [this preprint publication](https://psyarxiv.com/63a4y/), which has been peer reviewed and accepted for publication in [Scientific Data](http://nature.com/sdata/). The full [BIDS specification](https://bids-specification.readthedocs.io/en/stable/) including MEG, iEEG, MRI and behavioral data is maintained and available online.
%
%% # Converting to BIDS - copying the data
%
% If your EEG data is in a format supported by BIDS, you only have to rename the files and organize them in the right directory structure. After that, you have to add the required metadata.
%
% Rather than starting with 10 real EEG recordings, we can download a publicly available EEG file and copy it multiple times, pretending that this file contains the EEG data for each of our subjects. You can download the dataset from the EEGLAB website at <https://sccn.ucsd.edu/mediawiki/images/9/9c/Eeglab_data.set>.
%
sub = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10'};

% for subject 3 the age is unknown, for subject 2 the sex is not specified
age = [11  96  nan 77  82  87  18 40  26  80];
sex = {'f' [] 'f' 'f' 'f' 'm' 'm' 'm' 'm' 'm'};

for subindx=1:numel(sub)

  cfg = [];
  cfg.method    = 'copy';
  cfg.datatype  = 'eeg';

  % specify the input file name, here we are using the same file for every subject
  cfg.dataset   = dccnpath('/home/common/matlab/fieldtrip/data/ftp/example/bids_eeg/original/Eeglab_data.set');

  % specify the output directory
  cfg.bidsroot  = fullfile(tempdir,'bids');
  cfg.sub       = sub{subindx};

  % specify the information for the participants.tsv file
  % this is optional, you can also pass other pieces of info
  cfg.participants.age = age(subindx);
  cfg.participants.sex = sex{subindx};

  % specify the information for the scans.tsv file
  % this is optional, you can also pass other pieces of info
  cfg.scans.acq_time = datestr(now, 'yyyy-mm-ddThh:MM:SS'); % according to RFC3339

  % specify some general information that will be added to the eeg.json file
  cfg.InstitutionName             = 'University of California San Diego';
  cfg.InstitutionalDepartmentName = 'Schwartz Center for Computational Neuroscience';
  cfg.InstitutionAddress          = '9500 Gilman Drive # 0559; La Jolla CA 92093, USA';

  % provide the mnemonic and long description of the task
  cfg.TaskName        = 'changedetection';
  cfg.TaskDescription = 'Subjects were responding as fast as possible upon a change in a visually presented stimulus.';

  % these are EEG specific
  cfg.eeg.PowerLineFrequency = 60;   % since recorded in the USA
  cfg.eeg.EEGReference       = 'M1'; % actually I do not know, but let's assume it was left mastoid

  data2bids(cfg);

end
%
%% # Converting to BIDS - converting the data
%
% If your EEG data is in a format that is not supported by BIDS, you can convert it to BrainVision format. Contrary to many proprietary EEG formats, the BrainVision format is a simple file format that is well documented and widely supported in multiple programming environments (MATLAB, Python, C/C++) and all sorts of software (free and commercial). This makes it a very good format for storing data such that it remains accessible in the future.
%
% For this you would use the same script as above, but rather than
%
cfg.method = 'copy';

% you would specify
%
cfg.method = 'convert';

% The BrainVision format consists of a separate header, data and marker file. This means that each dataset now consists of three files, plus the TSV and JSON sidecar files for BIDS.
%
%% # Multiple sessions and/or runs
%
% With BIDS you can also organize data that was recorded in multiple sessions, i.e. lab visits, for example before and after a treatment or a night sleep. You can also specify multiple runs, i.e. subsequent recordings during the same session. Converting such data would go like this
%
sub = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10'};
ses = {'before', 'after'};
run = [1 2 3]; % this must be numeric

for subindx=1:numel(sub)
  for sesindx=1:numel(ses)
    for runindx=1:numel(run)

      cfg = [];
      cfg.method    = 'convert';
      cfg.datatype  = 'eeg';

      % specify the input file name, here we are using the same file for every subject
      cfg.dataset   = dccnpath('/home/common/matlab/fieldtrip/data/ftp/example/bids_eeg/original/Eeglab_data.set');

      % specify the output directory
      cfg.bidsroot  = fullfile(tempdir,'bids');
      cfg.sub       = sub{subindx};
      cfg.ses       = ses{sesindx};
      cfg.run       = run(runindx);

      data2bids(cfg);

    end % for run
  end % for ses
end % for sub
%
% Information that is session specific would be specified in `cfg.scans` and go in the `scans.tsv` file. Information that varies from one run to the next run (such as the sedation level, or the behavioral performance during the run) can also be added to the `scans.tsv` file.
%
% Information that is subject specific and identical for all recordings goes via `cfg.participants` into the `participants.tsv` file.
%
%% # Finalize the BIDS dataset
%
% There are some things which are not so conveniently implemented as a MATLAB script, for example filling out all details in the top-level `dataset_description.json` file, adding a |README| file and updating the |CHANGES| file. These are things you would usually do using a regular text editor after converting all data. You should also use the [bids-validator](http://github.com/bids-standard/bids-validator/) to check compliance. Note that the online validator does not require you to upload any data; it runs as JavaScript in your local browser, so the files stay on your computer.
%
% In converting data from one to another format there is always a chance of loosing information. Therefore you should store the data files in their original format in the `bids/sourcedata` directory. This also applies to auxiliary files, lab notes, presentation log files and other information that you do not want to loose.
%
% Furthermore, it is recommended to store the script that you used for converting the data in the `bids/code`. See the [BIDS specification](https://bids-specification.readthedocs.io/en/stable/) for further details.
%
%% # Concluding remarks
%
% In this example it all looks very simple, which is because the data files are perfectly identical and there are no exceptions. Usually challenge arises due to inconsistencies between recordings or subjects. Although these inconsistencies are not be part of the experimental protocol, stuff happens. For example
%
%* one subject was recorded with different acquisition settings
%* one subject was recorded with another EEG cap
%* for one subject some channels were swapped around
%* for one subject the experiment was paused, resulting in two recording files
%* for one subject the triggers were recorded incorrectly
%* etc.
%
% As a rule of thumb - if you have few exceptions, better don't try to make the scripts above too complex, but deal with them manually.
%
% If you have many exceptions of a similar type, it is worthwhile to invest into making these scripts smarter to automate the exception handling.
%
% In reusing the data, either by yourself, your (future) colleagues in your lab, or people outside your lab, this type of information is very relevant. Although it can be frustrating to deal with the lack of proper documentation and the inconsistencies in your data when converting to BIDS, it actually reveals that these aspects need to be represented and documented properly!
