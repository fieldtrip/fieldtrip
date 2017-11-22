function datasets = ref_datasets

% We have collected a whole bunch of datasets in the original acquisition
% systems file format that can be used for testing whether FT functions
% work properly for other types of data and not only for MEG data
% that was recorded on our ctf275 system.  These test datasetss (files,
% directories) are located on a shared network disk and are organized
% in a structure according to PATH/TYPE/DATATYPE/datafile or directory
% 
% where
%   PATH = somewhere on home/common
%   TYPE = eeg, meg, lfp
%   DATATYPE = for example ctf64, ctf151, ctf275

origdir = dccnpath('/home/common/matlab/fieldtrip/data/test/');

datasets = [];

% eeg dataformats
datasets(end+1).origdir = origdir; 
datasets(end).type      = 'eeg';
datasets(end).datatype  = 'bdf';
datasets(end).filename  = 'Newtest17-256.bdf';
datasets(end).numeeg    = nan;
datasets(end).nummeg    = nan;
datasets(end).numeog    = nan;
datasets(end).numecg    = nan;
datasets(end).numemg    = nan;
datasets(end).senstype  = 'unknown';

datasets(end+1).origdir = origdir; 
datasets(end).type      = 'eeg';
datasets(end).datatype  = 'brainvision';
datasets(end).filename  = 'Mischa.vhdr';
datasets(end).numeeg    = nan;
datasets(end).nummeg    = nan;
datasets(end).numeog    = 2;
datasets(end).numecg    = 2;
datasets(end).numemg    = 2;
datasets(end).senstype  = 'unknown';

datasets(end+1).origdir = origdir; 
datasets(end).type      = 'eeg';
datasets(end).datatype  = 'edf';
datasets(end).filename  = '0601_s.edf';
datasets(end).numeeg    = nan;
datasets(end).nummeg    = nan;
datasets(end).numeog    = nan;
datasets(end).numecg    = nan;
datasets(end).numemg    = nan;
datasets(end).senstype  = 'eeg';

datasets(end+1).origdir = origdir; 
datasets(end).type      = 'eeg';
datasets(end).datatype  = 'neuroscan16';
datasets(end).filename  = '0500.cnt';
datasets(end).numeeg    = 58;
datasets(end).nummeg    = nan;
datasets(end).numeog    = nan;
datasets(end).numecg    = nan;
datasets(end).numemg    = nan;
datasets(end).senstype  = 'ext1020';% type = ft_senstype(data.label)

datasets(end+1).origdir = origdir; 
datasets(end).type      = 'eeg';
datasets(end).datatype  = 'neuroscan32';
datasets(end).filename  = 'Subject1_MP.cnt';
datasets(end).numeeg    = 62;
datasets(end).nummeg    = nan;
datasets(end).numeog    = nan;
datasets(end).numecg    = nan;
datasets(end).numemg    = nan;
datasets(end).senstype  = 'ext1020';% type = ft_senstype(data.label)

% meg dataformats
datasets(end+1).origdir = origdir; 
datasets(end).type      = 'meg';
datasets(end).datatype  = 'bti148';
datasets(end).filename  = 'c,rfhp0.1Hz.m4d';
datasets(end).numeeg    = nan;
datasets(end).nummeg    = 148;
datasets(end).numeog    = nan;
datasets(end).numecg    = nan;
datasets(end).numemg    = nan;
datasets(end).senstype  = 'bti148';

datasets(end+1).origdir = origdir; 
datasets(end).type      = 'meg';
datasets(end).datatype  = 'bti248';
datasets(end).filename  = 'e,rfDC';
datasets(end).numeeg    = nan;
datasets(end).nummeg    = 248;
datasets(end).numeog    = nan;
datasets(end).numecg    = nan;
datasets(end).numemg    = nan;
datasets(end).senstype  = 'bti248';

datasets(end+1).origdir = origdir;
datasets(end).type      = 'meg';
datasets(end).datatype  = 'bti248grad';
datasets(end).filename  = 'e,rfhp1.0Hz,COH';
datasets(end).numeeg    = nan;
datasets(end).nummeg    = 248;
datasets(end).numeog    = nan;
datasets(end).numecg    = nan;
datasets(end).numemg    = nan;
datasets(end).senstype  = 'bti248grad';% type = ft_senstype(data.label)

datasets(end+1).origdir = origdir; 
datasets(end).type      = 'meg';
datasets(end).datatype  = 'ctf151';
datasets(end).filename  = 'Subject01.ds';
datasets(end).numeeg    = nan;
datasets(end).nummeg    = 151;
datasets(end).numeog    = 1;
datasets(end).numecg    = nan;
datasets(end).numemg    = nan;
datasets(end).senstype  = 'ctf151';

% % Temporary test for bug1976. Should not normally be executed
% datasets(end+1).origdir = origdir; 
% datasets(end).type      = 'meg';
% datasets(end).datatype  = 'ctf151_1trial1time';
% datasets(end).filename  = 'Subject01.ds';
% datasets(end).numeeg    = nan;
% datasets(end).nummeg    = 151;
% datasets(end).numeog    = 1;
% datasets(end).numecg    = nan;
% datasets(end).numemg    = nan;
% datasets(end).senstype  = 'ctf151';

datasets(end+1).origdir = origdir; 
datasets(end).type      = 'meg';
datasets(end).datatype  = 'ctf275';
datasets(end).filename  = 'A0132_Aud-Obj-Recognition_20051115_02.ds';
datasets(end).numeeg    = nan;
datasets(end).nummeg    = 275;
datasets(end).numeog    = nan;
datasets(end).numecg    = nan;
datasets(end).numemg    = nan;
datasets(end).senstype  = 'ctf275';

datasets(end+1).origdir = origdir; 
datasets(end).type      = 'meg';
datasets(end).datatype  = 'ctf64';
datasets(end).filename  = 'Wat123r1raw.ds';
datasets(end).dataformat = 'ctf_old';
datasets(end).numeeg    = nan;
datasets(end).nummeg    = 64;
datasets(end).numeog    = nan;
datasets(end).numecg    = nan;
datasets(end).numemg    = nan;
datasets(end).senstype  = 'ctf64';

datasets(end).numecg    = nan;
datasets(end+1).origdir = origdir; 
datasets(end).type      = 'meg';
datasets(end).datatype  = 'itab153';
datasets(end).filename  = 'srgcst85_0105.raw';
datasets(end).numeeg    = nan;
datasets(end).nummeg    = 153;
datasets(end).numeog    = nan;
datasets(end).numecg    = nan;
datasets(end).numemg    = nan;
datasets(end).senstype  = 'itab153';

datasets(end+1).origdir = origdir; 
datasets(end).type      = 'meg';
datasets(end).datatype  = 'itab28';
datasets(end).filename  = 'gibb0101.raw';
datasets(end).numeeg    = nan;
datasets(end).nummeg    = nan;
datasets(end).numeog    = nan;
datasets(end).numecg    = nan;
datasets(end).numemg    = nan;
datasets(end).senstype  = 'itab28';

datasets(end+1).origdir = origdir; 
datasets(end).type      = 'meg';
datasets(end).datatype  = 'itab28_old';
datasets(end).filename  = 'gibb0101.raw';
datasets(end).numeeg    = nan;
datasets(end).nummeg    = nan;
datasets(end).numeog    = nan;
datasets(end).numecg    = nan;
datasets(end).numemg    = nan;
datasets(end).senstype  = 'itab28';

datasets(end+1).origdir = origdir; 
datasets(end).type      = 'meg';
datasets(end).datatype  = 'neuromag122';
datasets(end).filename  = 'jg_single_01raw.fif';
datasets(end).numeeg    = nan;
datasets(end).nummeg    = 122;
datasets(end).numeog    = nan;
datasets(end).numecg    = nan;
datasets(end).numemg    = nan;
datasets(end).senstype  = 'neuromag122';

datasets(end+1).origdir = origdir; 
datasets(end).type      = 'meg';
datasets(end).datatype  = 'neuromag306';
datasets(end).filename  = 'raw.fif';
datasets(end).numeeg    = 2;
datasets(end).nummeg    = 306;
datasets(end).nummegmag = 102;
datasets(end).nummeggrad = 204;
datasets(end).numeog    = nan;
datasets(end).numecg    = nan;
datasets(end).numemg    = nan;
datasets(end).senstype  = 'neuromag306';

datasets(end+1).origdir = origdir; 
datasets(end).type      = 'meg';
datasets(end).datatype  = 'yokogawa160';
datasets(end).filename  = 'Continuous1.con';
datasets(end).numeeg    = nan;
datasets(end).nummeg    = 160;
datasets(end).nummegmag = nan;
datasets(end).nummeggrad = 160;
datasets(end).numeog    = nan;
datasets(end).numecg    = nan;
datasets(end).numemg    = nan;
datasets(end).senstype  = 'yokogawa160';

% lfp dataformats 
datasets(end+1).origdir = origdir; 
datasets(end).type      = 'lfp';
datasets(end).datatype  = 'neuralynx';
datasets(end).filename  = '2005-2-10_11-57-40';
datasets(end).numeeg    = nan;
datasets(end).nummeg    = nan;
datasets(end).numeog    = nan;
datasets(end).numecg    = nan;
datasets(end).numemg    = nan;
datasets(end).senstype  = 'unknown';

% first trial has lot's of NaNs in raw data
% datasets(end+1).origdir = origdir; 
% datasets(end).type     = 'lfp';
% datasets(end).datatype = 'plexon';
% datasets(end).filename = 'p213parall.nex';
