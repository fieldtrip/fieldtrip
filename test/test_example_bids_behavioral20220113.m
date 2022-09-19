function test_example_bids_behavioral

% MEM 4gb
% WALLTIME 00:10:00

%
%% Converting an example behavioral dataset for sharing in BIDS
%
% Although [BIDS](https://bids.neuroimaging.io) was initiated for brain imaging data, possibly recorded while the participant was performing a task including stimuli and responses, it also allows for representing behavioural data that is recorded on its own. For example, you might run an experiment with NBS Presentation, PsychToolbox, or PsychoPy, in which you are only recording the behavioural responses to the stimuli.
%
% The [BIDS specification](https://bids-specification.readthedocs.io/) explains how behavioral events (and their timing) can be represented.
%
%% # The data to consider
%
% The raw data that is acquired in a typical behavioural experiment consists of one or multiple ASCII log files on disk. Additional data that is not stored in a systematic way may include demographic information of the participant, the details of the task and the instruction, but also the detailed stimulus material (e.g., wav and png files). Some results might already be derived during the experiment, such as whether responses were correct or incorrect (coded in the ASCII file), but might also be processed offline (e.g., rejecting trials with a response that is too fast or slow, excluding trials at the start of the experiment, or immediately following an incorrect trial). This processing is limited to some a priori design choices, and excludes the processing and analysis related to the research question.
%
% The data that you want to store for analysis and/or reuse combines the original raw ASCII log files, the non-systematic subject and task information, and potentially also the results of minimal processing. The example below shows how MATLAB and **[data2bids](https://github.com/fieldtrip/fieldtrip/blob/release/data2bids.m)** are used to organize the data in a reusable format.
%
%% # Examples
%
% When using NBS Presentation, you can write the timing of stimuli and responses in two different ways: using the standard Presentation .log file, or using your own code that is part of your presentation script. When using PsychToolbox or PsychoPy, but also using Presentation, you can write your events to a custom tabular format. Both the Presentation log file and a custom tabular format are demonstrated below.
%
% The example consists of data from 2 subjects, each participating in 3 experiments that were part of a single session. In one experiment the Presentation log file was stored, in the other two events were stored in a custom log file.
%
% All data for the following examples is available from our [FTP server](ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/example/bids_presentation/).
%
% The examples includes the original log files under the “original” directory. It also has a copy of the script to do the conversion under “code”. The reorganized data is under the “bids” directory. According to [the documentation](https://bids-specification.readthedocs.io/en/stable/02-common-principles.html#source-vs-raw-vs-derived-data) the original data can be added to the BIDS dataset under the “sourcedata” directory and code can be added to the “code” directory. This way no information is lost and the conversion/reorganization is fully reproducible.
%
%%

% The original log files contain some information:
% - the first part of the name contains the subject identifier: P31237 and P23671
% - the second part describes the task
% - the extension shows whether it is the original NBS Presentation log file itself (.log) or the custom text file written by the Presentation script (.txt)

%sourcepath = './original/sourcedata';
sourcepath = dccnpath('/home/common/matlab/fieldtrip/data/ftp/example/bids_presentation/original');
%targetpath = './bids';
targetpath = fullfile(tempdir, 'bids');

%%
% general information that applies to all

general= [];
general.InstitutionName                 = 'Radboud University';
general.InstitutionalDepartmentName     = 'Donders Institute for Brain, Cognition and Behaviour';
general.InstitutionAddress              = 'Kapittelweg 29, 6525 EN, Nijmegen, The Netherlands';
% this goes into dataset_description.json
general.dataset_description.Name        = 'Presentation example';
general.dataset_description.BIDSVersion = 'unofficial extension';
% specify the type and target location
general.bidsroot = targetpath;
general.datatype = 'events';
%
%% ## Presentation log files
%
% the first two files are in Presentation .log format

filename = {
  'P31237_task1-MotorTaskEv_left.log'
  'P23671_task1-MotorTaskEv_right.log'
  };

subject = {'P31237', 'P23671'};
task = {'motor', 'motor'};
age = [45, 52];
handedness = {'right', 'right'};
clinical_state = {'parkinson', 'control'};

for i=1:numel(subject)

  cfg = general; % copy the general fields
  cfg.sub = subject{i};
  cfg.task = task{i};

  % you can use this to add information to the participants.tsv table
  cfg.participants.age = age(i);
  cfg.participants.handedness = handedness(i);
  cfg.participants.clinical_state = clinical_state{i};

  cfg.dataset = fullfile(sourcepath, filename{i});

  % instead of specifying cfg.dataset, you can also use FT_READ_EVENT and pass the events like this
  % cfg.events = ft_read_event(fullfile(sourcepath, filename{i}));

  data2bids(cfg);
end
%
%% ## Custom log files
%
% The standard MATLAB [readtable](https://www.mathworks.com/help/matlab/ref/readtable.html) function is very useful to read tabular data. It automatically figures out whether there are header lines to skip, what the column headings are (if any), and what the separator between the columns is. This function returns a MATLAB table object, which is similar to a 2D cell-array, but it includes column headings.
%
% BIDS requires that the events include the onset and the duration (both in seconds) as the first and second column in the events.tsv file. Here we compute the onset and duration using the fixation and response time, which are two columns in the log file.
%
% the subsequent files are in a custom tabular format

filename = {
  'P31237_prac1_logfile.txt'
  'P23671_prac1_logfile.txt'
  'P31237_task1_logfile.txt'
  'P23671_task1_logfile.txt'
  };

subject = {'P31237', 'P23671', 'P31237', 'P23671'};
task = {'practice', 'practice', 'unknown', 'unknown'};

for i=1:numel(filename)

  cfg = general; % copy the general fields
  cfg.sub = subject{i};
  cfg.task = task{i};

  % read the ascii log file
  log = readtable(fullfile(sourcepath, filename{i}));

  % add the onset and duration (both in seconds)
  % the Presentation software uses time stamps of 0.1 milliseconds
  log.onset = (log.Fixation_Time)/1e4;
  log.duration = (log.Response_Time - log.Fixation_Time)/1e4;
  log.duration(log.Response_Time==0) = nan;

  cfg.events = log;
  data2bids(cfg);
end
