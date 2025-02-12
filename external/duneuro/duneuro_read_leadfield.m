function cfg = duneuro_read_leadfield(cfg)
% bst_read_binary_leadfield_matrix(filename)
% read the binary output from the duneuro application.
% could be used to read either the eeg/meg transfer matrix or the eeg/meg
% leadfield mtrix.
% usage :

% binary outpuf from duneuro :
% 'eeg_transfer.dat'
% 'meg_transfer.dat'
% 'meg_lf.dat'
% 'eeg_lf.dat'

switch cfg.modality
  case 'eeg'
    filename = fullfile(cfg.outputpath, 'eeg_lf.dat');
    data = duneuro_read_binary(filename);
    cfg.eeg.lf = data';
  case 'meg'
    filename = fullfile(cfg.outputpath, 'meg_lf.dat');
    data = duneuro_read_binary(filename);
    cfg.meg.Bs = data';
  case 'meeg'
    filename = fullfile(cfg.outputpath, 'meg_lf.dat');
    data = duneuro_read_binary(filename);
    cfg.meg.Bs = data';
    filename = fullfile(cfg.outputpath, 'eeg_lf.dat');
    data = duneuro_read_binary(filename);
    cfg.eeg.lf = data';
end

% to be completed with other modalities seeg, ieeg ...

