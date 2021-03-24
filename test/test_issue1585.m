function test_issue1585

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_read_header ft_read_data ft_read_event

%%

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1585'))

% cfg = [];
% cfg.dataset      = 'eeglab_data.mat';
% cfg.datatype     = 'eeg';

% cfg.headerformat = 'eeglab_set';
% cfg.dataformat   = 'eeglab_set';
% cfg.eventformat  = 'eeglab_set';
%
% cfg.method    = 'convert';
% cfg.bidsroot  = './bids';
% cfg.sub       = '001';
% cfg.task      = 'whatever';
% data2bids(cfg);

% after creating the BIDS dataset, I used a text editor to remove a row from the channels.tsv

%%

try
  % this should give an error
  hdr = ft_read_header('bids/sub-001/eeg/sub-001_task-whatever_eeg.eeg', 'readbids', 'yes');
  dat = ft_read_data('bids/sub-001/eeg/sub-001_task-whatever_eeg.eeg', 'readbids', 'yes');
  evt = ft_read_event('bids/sub-001/eeg/sub-001_task-whatever_eeg.eeg', 'readbids', 'yes');
  failed = false;
catch
  failed = true;
end
assert(failed, 'it should have failed');

% this should give a warning
hdr = ft_read_header('bids/sub-001/eeg/sub-001_task-whatever_eeg.eeg', 'readbids', 'ifmakessense');
dat = ft_read_data('bids/sub-001/eeg/sub-001_task-whatever_eeg.eeg', 'readbids', 'ifmakessense');
evt = ft_read_event('bids/sub-001/eeg/sub-001_task-whatever_eeg.eeg', 'readbids', 'ifmakessense');

% this should work without any complaints
hdr = ft_read_header('bids/sub-001/eeg/sub-001_task-whatever_eeg.eeg', 'readbids', 'no');
dat = ft_read_data('bids/sub-001/eeg/sub-001_task-whatever_eeg.eeg', 'readbids', 'no');
evt = ft_read_event('bids/sub-001/eeg/sub-001_task-whatever_eeg.eeg', 'readbids', 'no');
