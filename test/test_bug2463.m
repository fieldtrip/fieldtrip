function test_bug2463

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_bug24673
% TEST ft_channelselection ft_chantype

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2463.mat'));

% data.label'
% ans =
%     'Fp1'    'Fp2'    'F3'     'F4'     'C3'     'C4'     'P3'     'P4'     'O1'     'O2'    'F7'    'F8'    'T7'    'T8'    'P7'    'P8'    'Fz'    'Cz'    'Pz'
%     'Oz'     'FC1'    'FC2'    'CP1'    'CP2'    'FC5'    'FC6'    'CP5'    'CP6'    'F1'    'F2'    'C1'    'C2'    'P1'    'P2'    'AF3'   'AF4'
%     'FC3'    'FC4'    'CP3'    'CP4'    'PO3'    'PO4'    'F5'     'F6'     'C5'     'C6'    'P5'    'P6'    'AF7'   'AF8'   'FT7'   'FT8'   'TP7'
%     'TP8'    'PO7'    'PO8'    'Fpz'    'AFz'    'CPz'    'POz'    'FCz'

% the problem is that only 21 from the 61 channels get selected, although the
% channels are clearly all EEG
assert(length(ft_channelselection('EEG', data.label))==61, 'not all 61 EEG channels were selected');

% the problem is also visible in ft_chantype
%
% cat(2, data.label, ft_chantype(data.label))
% ans = 
%     'Fp1'    'eeg'    
%     'Fp2'    'eeg'    
%     'F3'     'eeg'    
%     'F4'     'eeg'    
%     'C3'     'eeg'    
%     'C4'     'eeg'    
%     'P3'     'eeg'    
%     'P4'     'eeg'    
%     'O1'     'eeg'    
%     'O2'     'eeg'    
%     'F7'     'eeg'    
%     'F8'     'eeg'    
%     'T7'     'eeg'    
%     'T8'     'eeg'    
%     'P7'     'eeg'    
%     'P8'     'eeg'    
%     'Fz'     'eeg'    
%     'Cz'     'eeg'    
%     'Pz'     'eeg'    
%     'Oz'     'eeg'    
%     'FC1'    'unknown'
%     'FC2'    'unknown'
%     'CP1'    'unknown'
%     'CP2'    'unknown'
%     'FC5'    'unknown'
%     'FC6'    'unknown'
%     'CP5'    'unknown'
%     'CP6'    'unknown'
%     'F1'     'unknown'
%     'F2'     'unknown'
%     'C1'     'unknown'
%     'C2'     'unknown'
%     'P1'     'unknown'
%     'P2'     'unknown'
%     'AF3'    'unknown'
%     'AF4'    'unknown'
%     'FC3'    'unknown'
%     'FC4'    'unknown'
%     'CP3'    'unknown'
%     'CP4'    'unknown'
%     'PO3'    'unknown'
%     'PO4'    'unknown'
%     'F5'     'unknown'
%     'F6'     'unknown'
%     'C5'     'unknown'
%     'C6'     'unknown'
%     'P5'     'unknown'
%     'P6'     'unknown'
%     'AF7'    'unknown'
%     'AF8'    'unknown'
%     'FT7'    'unknown'
%     'FT8'    'unknown'
%     'TP7'    'unknown'
%     'TP8'    'unknown'
%     'PO7'    'unknown'
%     'PO8'    'unknown'
%     'Fpz'    'eeg'    
%     'AFz'    'unknown'
%     'CPz'    'unknown'
%     'POz'    'unknown'
%     'FCz'    'unknown'


