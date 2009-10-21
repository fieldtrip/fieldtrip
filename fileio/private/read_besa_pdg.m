function [bad_flags, epochs, rejected, trial_pos_ms, trial_cond_ind, scanned_range] = read_besa_pdg(filename);

% Reads BESA *.pdg file and finds out what trials were rejected
% filename - PDG file name
% rejected - array of flags: 1 if the trial was rejected, 0 otherwise
% bad_flags - flags marking bad channels. 0 - channel was good 
% trial_pos_ms - porition of trigger in ms for each trial
% trial_cond_ind - index of condition to which each trial belongs
% epochs - matrix describing the epochs relative to trigger
% scanned_range - start and end (in ms) of the data range scanned

% Copyright (C) 2005, Vladimir Litvak 6/4/05
%
% $Log: read_besa_pdg.m,v $
% Revision 1.1  2009/01/14 09:24:45  roboos
% moved even more files from fileio to fileio/privtae, see previous log entry
%
% Revision 1.2  2005/07/28 14:08:06  roboos
% converted from partially dos/unix into pure unix (i.e. removed CRs)
%
% Revision 1.1  2005/07/28 13:59:09  roboos
% renamed read_pdg.m into read_besa_pdg.m
%
% Revision 1.2  2005/07/28 13:58:10  roboos
% added copyrights, added log
% updated to latest version from Vladimir (1/7/05)
%

data_type=1; % needs to be 1 for EEG, 2 for MEG

pdg_file=fopen(filename, 'r');

epochs=[];
% Looks for 'Epochs' part in the file
while (~feof(pdg_file) & isempty(findstr(fgetl(pdg_file), '[Epochs]'))) end;
if feof(pdg_file)
    error('Illegal format');
end

% Reads epochs information
curr_line=fgetl(pdg_file);
while ~all(isspace(curr_line))
    epochs=[epochs ; str2num(curr_line)];
    curr_line=fgetl(pdg_file);
end

% Looks for 'Thresholds' part in the file
while (~feof(pdg_file) & isempty(findstr(fgetl(pdg_file), '[Thresholds]'))) end;
if feof(pdg_file)
    error('Illegal format');
end

% Skips the part having to do with GUI color scaling which is irrelevant
while (~feof(pdg_file) & isempty(findstr(fgetl(pdg_file), 'AUTO_REJECT'))) end;
if feof(pdg_file)
    error('No auto rejection was performed');
end

% Reads thresholds information
thresholds=[];
for i=1:3
    thresholds=[thresholds ; str2num(fgetl(pdg_file))];
end

% Looks for 'ArtifactScan' part in the file
while (~feof(pdg_file) & isempty(findstr(fgetl(pdg_file), '[ArtifactScan]'))) end;
if feof(pdg_file)
    error('Illegal format');
end

% Not very useful
chan_numbers=str2num(fgetl(pdg_file));

% Marks rejected channels. Non-zero entries mean that channel was rejected
bad_flags=str2num(fgetl(pdg_file));

% If only bad flags are necessary there is no need to read the long part with all the trials
if nargout<3
    fclose(pdg_file);
    return;
end

% Should used for finding out to which condition each trial belongs
cond_num_ind=str2num(fgetl(pdg_file));

scanned_range=str2num(fgetl(pdg_file))./1000;

% Marks the condition of each trial
trial_cond_ind=[];

% Trigger position in us
trial_pos=[];

% Flag marking if the trial was rejected.
rejected=[];
while ~feof(pdg_file)
    tmp=str2num(fgetl(pdg_file));
    trial_cond_ind=[trial_cond_ind tmp(1)];
    trial_pos=[trial_pos tmp(2)];
    
    amp_trial=str2num(fgetl(pdg_file));
    var_trial=str2num(fgetl(pdg_file));
    grad_trial=str2num(fgetl(pdg_file));
    
    rejected=[rejected any([...
                any(amp_trial(~bad_flags)>thresholds(1,data_type)) ... 
                any(var_trial(~bad_flags)<thresholds(2,data_type))...
                any(grad_trial(~bad_flags)>thresholds(3,data_type))])];
    
end

% Translates position to ms.
trial_pos_ms=trial_pos/1000;

fclose(pdg_file);
