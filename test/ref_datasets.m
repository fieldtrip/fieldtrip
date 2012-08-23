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

if isunix
  origdir = '/home/common/matlab/fieldtrip/data/test/';
elseif ispc
  origdir = fullfile('H:', 'common', 'matlab', 'fieldtrip', 'data', 'test');
end

% eeg dataformats
datasets(1).origdir  = origdir; 
datasets(1).type     = 'eeg';
datasets(1).datatype = 'bdf';
datasets(1).filename = 'Newtest17-256.bdf';

datasets(end+1).origdir  = origdir; 
datasets(end).type     = 'eeg';
datasets(end).datatype = 'brainvision';
datasets(end).filename = 'Mischa.vhdr';

datasets(end+1).origdir  = origdir; 
datasets(end).type     = 'eeg';
datasets(end).datatype = 'edf';
datasets(end).filename = '0601_s.edf';

datasets(end+1).origdir  = origdir; 
datasets(end).type     = 'eeg';
datasets(end).datatype = 'neuroscan16';
datasets(end).filename = '0500.cnt';

datasets(end+1).origdir  = origdir; 
datasets(end).type     = 'eeg';
datasets(end).datatype = 'neuroscan32';
datasets(end).filename = 'Subject1_MP.cnt';

% meg dataformats
datasets(end+1).origdir  = origdir; 
datasets(end).type     = 'meg';
datasets(end).datatype = 'bti148';
datasets(end).filename = 'c,rfhp0.1Hz.m4d';

datasets(end+1).origdir  = origdir; 
datasets(end).type     = 'meg';
datasets(end).datatype = 'bti248';
datasets(end).filename = 'e,rfDC';

datasets(end+1).origdir  = origdir;
datasets(end).type     = 'meg';
datasets(end).datatype = 'bti248grad';
datasets(end).filename = 'e,rfhp1.0Hz,COH';

datasets(end+1).origdir  = origdir; 
datasets(end).type     = 'meg';
datasets(end).datatype = 'ctf151';
datasets(end).filename = 'Subject01.ds';

datasets(end+1).origdir  = origdir; 
datasets(end).type     = 'meg';
datasets(end).datatype = 'ctf275';
datasets(end).filename = 'A0132_Aud-Obj-Recognition_20051115_02.ds';

datasets(end+1).origdir  = origdir; 
datasets(end).type     = 'meg';
datasets(end).datatype = 'ctf64';
datasets(end).filename = 'Wat123r1raw.ds';
datasets(end).dataformat = 'ctf_old';

datasets(end+1).origdir  = origdir; 
datasets(end).type     = 'meg';
datasets(end).datatype = 'itab153';
datasets(end).filename = 'srgcst85_0105.raw';

datasets(end+1).origdir  = origdir; 
datasets(end).type     = 'meg';
datasets(end).datatype = 'itab28';
datasets(end).filename = 'gibb0101.raw';

datasets(end+1).origdir  = origdir; 
datasets(end).type     = 'meg';
datasets(end).datatype = 'itab28_old';
datasets(end).filename = 'gibb0101.raw';

datasets(end+1).origdir  = origdir; 
datasets(end).type     = 'meg';
datasets(end).datatype = 'neuromag122';
datasets(end).filename = 'jg_single_01raw.fif';

datasets(end+1).origdir  = origdir; 
datasets(end).type     = 'meg';
datasets(end).datatype = 'neuromag306';
datasets(end).filename = 'raw.fif';

datasets(end+1).origdir  = origdir; 
datasets(end).type     = 'meg';
datasets(end).datatype = 'yokogawa160';
datasets(end).filename = 'Continuous1.con';

% lfp dataformats
% datasets(end+1).origdir  = origdir; 
% datasets(end).type     = 'lfp';
% datasets(end).datatype = 'neuralynx';
% datasets(end).filename = '2005-2-10_12-11-19';

datasets(end+1).origdir  = origdir; 
datasets(end).type     = 'lfp';
datasets(end).datatype = 'plexon';
datasets(end).filename = 'p213parall.nex';
