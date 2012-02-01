function dataset = datasets

origdir = '/home/common/matlab/fieldtrip/data/test/latest/';

% eeg dataformats
dataset(1).origdir  = origdir; 
dataset(1).type     = 'eeg/';
dataset(1).datatype = 'bdf';
dataset(1).filename = 'Newtest17-256.bdf';

dataset(end+1).origdir  = origdir; 
dataset(end).type     = 'eeg/';
dataset(end).datatype = 'brainvision';
dataset(end).filename = 'Mischa.vhdr';

dataset(end+1).origdir  = origdir; 
dataset(end).type     = 'eeg/';
dataset(end).datatype = 'edf';
dataset(end).filename = '0601_s.edf';

dataset(end+1).origdir  = origdir; 
dataset(end).type     = 'eeg/';
dataset(end).datatype = 'neuroscan16';
dataset(end).filename = '0500.cnt';

dataset(end+1).origdir  = origdir; 
dataset(end).type     = 'eeg/';
dataset(end).datatype = 'neuroscan32';
dataset(end).filename = 'Subject1_MP.cnt';

% meg dataformats
dataset(end+1).origdir  = origdir; 
dataset(end).type     = 'meg/';
dataset(end).datatype = 'bti148';
dataset(end).filename = 'c,rfhp0.1Hz.m4d';

dataset(end+1).origdir  = origdir; 
dataset(end).type     = 'meg/';
dataset(end).datatype = 'bti248';
dataset(end).filename = 'e,rfDC';

dataset(end+1).origdir  = origdir;
dataset(end).type     = 'meg/';
dataset(end).datatype = 'bti248grad';
dataset(end).filename = 'e,rfhp1.0Hz,COH';

dataset(end+1).origdir  = origdir; 
dataset(end).type     = 'meg/';
dataset(end).datatype = 'ctf151';
dataset(end).filename = 'Subject01.ds';

dataset(end+1).origdir  = origdir; 
dataset(end).type     = 'meg/';
dataset(end).datatype = 'ctf275';
dataset(end).filename = 'A0132_Aud-Obj-Recognition_20051115_02.ds';

dataset(end+1).origdir  = origdir; 
dataset(end).type     = 'meg/';
dataset(end).datatype = 'ctf64';
dataset(end).filename = 'Wat123r1raw.ds';
dataset(end).dataformat = 'ctf_old';

dataset(end+1).origdir  = origdir; 
dataset(end).type     = 'meg/';
dataset(end).datatype = 'itab153';
dataset(end).filename = 'srgcst85_0105.raw';

dataset(end+1).origdir  = origdir; 
dataset(end).type     = 'meg/';
dataset(end).datatype = 'itab28';
dataset(end).filename = 'gibb0101.raw';

dataset(end+1).origdir  = origdir; 
dataset(end).type     = 'meg/';
dataset(end).datatype = 'itab28_old';
dataset(end).filename = 'gibb0101.raw';

dataset(end+1).origdir  = origdir; 
dataset(end).type     = 'meg/';
dataset(end).datatype = 'neuromag122';
dataset(end).filename = 'jg_single_01raw.fif';

dataset(end+1).origdir  = origdir; 
dataset(end).type     = 'meg/';
dataset(end).datatype = 'neuromag306';
dataset(end).filename = 'raw.fif';

dataset(end+1).origdir  = origdir; 
dataset(end).type     = 'meg/';
dataset(end).datatype = 'yokogawa160';
dataset(end).filename = 'Continuous1.con';

% lfp dataformats
% dataset(end+1).origdir  = origdir; 
% dataset(end).type     = 'lfp/';
% dataset(end).datatype = 'neuralynx';
% dataset(end).filename = '2005-2-10_12-11-19';

dataset(end+1).origdir  = origdir; 
dataset(end).type     = 'lfp/';
dataset(end).datatype = 'plexon';
dataset(end).filename = 'p213parall.nex';
