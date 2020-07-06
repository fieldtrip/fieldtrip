% eeg_options() - eeglab option script 
%
% Note: DO NOT EDIT, instead use pop_editoptions() or the menu
%       /File/Maximize memory in EEGLAB gui

% STUDY and file options (set these checkboxes if you intend to work with studies) 
option_storedisk     = 0 ;  % If set, keep at most one dataset in memory. This allows processing hundreds of datasets within studies.
option_savetwofiles  = 1 ;  % If set, save not one but two files for each dataset (header and data). This allows faster data loading in studies.
option_saveversion6  = 1 ;  % If set, write Matlab files in Matlab v6.5 (max compatibility). If not, write files in Matlab v7.3 (larger than 2Gb).
% Memory options 
option_single        = 1 ;  % If set, use single precision under Matlab 7.x. This saves RAM but can lead to rare numerical imprecisions. 
option_memmapdata    = 0 ;  % If set, use memory mapped array under Matlab 7.x. This may slow down some computation (beta).
option_eegobject     = 0 ;  % If set, use the EEGLAB EEG object instead of the standard EEG structure (beta). 
% ICA options 
option_computeica    = 1 ;  % If set, precompute ICA activations. This requires more RAM but allows faster plotting of component activations. 
option_scaleicarms   = 1 ;  % If set, scale ICA component activities to RMS (Root Mean Square) in microvolt (recommended).
% Folder options
option_rememberfolder = 1 ;  % If set, when browsing to open a new dataset assume the folder/directory of previous dataset.
% Toolbox options
option_donotusetoolboxes = 0 ;  % If set, do not use Matlab additional toolboxes functions even if they are present (need to restart EEGLAB).
% EEGLAB connectivity and support
option_checkversion = 1 ;  % If set, check for new version of EEGLAB at startup.
option_chat = 0 ;  % If set, enable EEGLAB chat - currently UCSD only - restart EEGLAB after changing that option.
