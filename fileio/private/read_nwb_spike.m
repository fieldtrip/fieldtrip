function [spike] = read_nwb_spike(filename)
% READ_NWB_SPIKE reads spike timestamps and waveforms (if present-not currently supported) 
% from NWB files and converts them to fieldtrip spike data format. 
%
% INPUT: filename = (Path and) name of the .nwb file
%
% OUTPUT: spike = FieldTrip spike structure
%
% Notes: 
% This function was written during the NWB hackathon May 2020. It is
% based on example data in .nwb format schema version 2.0.1: 
% https://gui.dandiarchive.org/#/file-browser/folder/5e6eb2b776569eb93f451f8d
%
% NWB is a complicated data format and under active development. We
% recommend to use the latest stable release of MatNWB from the github
% page: https://github.com/NeurodataWithoutBorders/matnwb/releases
% and familiarize yourself with the use of generateCore():
% https://neurodatawithoutborders.github.io/matnwb
%
% With util.getSchemaVersion(file.nwb) the nwb file version can be
% querried. It may be necessary to replace the files in ..\matnwb\nwb-schema\core
% with the files from the nwb-schema version the file was created in from
% ..\nwb-schema\core. 
% Nwb-schemas can be obtained from here: 
% https://github.com/NeurodataWithoutBorders/nwb-schema/releases
% 
% -----
% Latest change: 01/06/2020


% Load data 
nwb = nwbRead(filename);
nwb_schema_v = util.getSchemaVersion(filename);
fprintf('\nReading NWB data from \n%s\n', filename);
fprintf('NWB schema version: %s\n', nwb_schema_v{1})



% Get spike sampling frequency, if not available get LFP sampling frequency
% instead and warn the user
Fs = 1; % if we cannot determine Fs, set it to 1 (timestamps will be equal to time)
foundin = 'not found in NWB file';
try
    % Check if Fs is saved with spike times
    if ismember({'resolution'}, fieldnames(nwb.units))
        Fs = nwb.units.resolution;
        foundin = 'units.resolution';
    % Check if Fs is saved with SpikeEventSeries 
    else
        processing_keys = nwb.processing.keys();
        if ismember({'ecephys'}, processing_keys)
            keys = nwb.processing.get('ecephys').nwbdatainterface.keys();
            for ikey=1:length(keys)
                time_series = nwb.processing.get('ecephys').nwbdatainterface.get(keys{ikey});
                if isa(time_series, 'types.core.SpikeEventSeries')
                    if ~isempty(time_series.starting_time_rate)
                        Fs = time_series.starting_time_rate;
                        foundin = 'SpikeEventSeries class';
                    end
                end
            end
            % if no spiking Fs exists, get LFP Fs
            if ismember({'LFP'}, keys) && strcmp(foundin, 'not found in NWB file')
                warning('Could not find spike Fs specifically. Looking for LFP Fs.');
                lfp_set = nwb.processing.get('ecephys').nwbdatainterface.get('LFP');
                fn = fieldnames(lfp_set);
                next_key = lfp_set.(fn{1}).keys();
                nwb_lfp = lfp_set.(fn{1}).get(next_key{1});
                Fs = nwb_lfp.starting_time_rate;
                foundin = 'LFP class';
            end
        end
    end
    if strcmp(foundin, 'not found in NWB file')
        warning('Could not find Fs in NWB file. Assuming Fs = 1')
    else
        fprintf('\nTrying to determine Fs from NWB file, setting to %i, found in: %s.', Fs, foundin);
    end
catch ME
    warning('Error while looking for sampling frequency from NWB file.');
    rethrow(ME)
end



% Get unit IDS, spike times and indeces from NWB file
%
nwb_unit_ids = nwb.units.id.data.load();
% spike times in s(?)
nwb_spike_times = nwb.units.spike_times.data.load();
% Index of start of unit in spike_times
nwb_spike_times_idx = nwb.units.spike_times_index.data.load();
Nunits = length(nwb_unit_ids);



% Check how long a recording is in hours (approximately)
dur_rec_in_h = (nwb_spike_times(nwb_spike_times_idx(1)) - nwb_spike_times(1)) / 3600;
fprintf('\nDuration of recording is %.2f hours\n', dur_rec_in_h)



% Trial info 
if ~isempty(nwb.intervals_trials)
    fprintf('\nFound trial information.\n');
    nwb_trial_ids = nwb.intervals_trials.id.data.load();
    % start and end time of trials
    nwb_trial_stime = nwb.intervals_trials.start_time.data.load();
    nwb_trial_etime = nwb.intervals_trials.stop_time.data.load();
    Ntrials = length(nwb_trial_ids);
else
    fprintf('\nNo trial information found\n.');
end



% Initialize FT fields for ft_datatype_spike
label = cell(1,Nunits);         % channel labels
timestamp = cell(1,Nunits);     % spike samples
time = cell(1,Nunits);          % spike times
trial = cell(1,Nunits);         % trial index (for each spike)
unit = cell(1,Nunits);          % unit index for each spike
trialtime = [];                 % time of start and end of trial

% Go through NWB units and fill in FT fields
fprintf('\nConverting spikes from NWB to FieldTrip:\n')
last_idx = 0;
for i = 1:Nunits

    if i == 1 || i == Nunits || mod(i, round(Nunits/10)) == 0
        fprintf('Processing unit %i \n', i)
    end
    
    % get current unit ID
    unit_id = nwb_unit_ids(i);
    
    % Map units to electrodes (add label)
    label{i} = ['unit', num2str(unit_id)];
    
    % Determine start and end of unit in nwb_spike_times
    start_idx = last_idx + 1;
    end_idx = nwb_spike_times_idx(i);
    last_idx = end_idx;
    
    % Add spike times for each unit
    time{i} = nwb_spike_times(start_idx:end_idx);
    
    % Add spike samples for each unit. 
    timestamp{i} = nwb_spike_times(start_idx:end_idx) * Fs;
    
    % Add physical measurement unit field - e.g. V (TO DO: Find unit in NWB file!)
    unit{i} = nan(size(time{i}));
    
    % Add trial structure if it exists
    if exist('nwb_trial_stime', 'var')
        trial{i} = nan(size(time{i}));
        for itrial = 1:Ntrials
            % find all spike times in trial itrial
            trl_idx = time{i}>=nwb_trial_stime(itrial) & time{i}<nwb_trial_etime(itrial);
            trial{i}(trl_idx) = itrial;
        end
        
        % Delete spikes that are not in trials? Otherwise filled with NaNs.
        % leave nans for now
                
        % make sure times are relative to event of interest in trial?
        % Unclear if that helps the user much - leave for them to handle
    end

end


% Add unit (i.e. cell/neuron) information and mappings to ft struct
nwb_unit_info = constrained_set_to_struct(nwb.units.vectordata, ...
                                                'from nwb.units.vectordata');

% Add stimulus info from NWB to ft struct
nwb_stim_info = constrained_set_to_struct(nwb.stimulus_presentation, ...
                                                'from nwb.stimulus_presentation');
                                            
% Add nwb trial information to ft struct
if exist('nwb_trial_stime', 'var')    
    nwb_trial_info = constrained_set_to_struct(nwb.intervals_trials.vectordata, ...
                                                'from nwb.intervals_trials.vectordata');
end

% Add nwb electrode information to ft struct
nwb_elec_info = constrained_set_to_struct(nwb.general_extracellular_ephys_electrodes.vectordata, ...
                                                'from nwb.general_extracellular_ephys_electrodes.vectordata');
elecs = nwb.general_extracellular_ephys_electrodes;
nwb_elec_info.elec_label = elecs.vectordata.get('group_name').data.load;
nwb_elec_info.elec_id = elecs.id.data.load;
if ~isempty(nwb.units.electrode_group)
    nwb_elec_info.unit_to_elec = {nwb.units.electrode_group.data.path};
end


% Create fieldtrip data structure
spike.label = label;
spike.timestamp = timestamp;
spike.time = time;
spike.unit = unit;
% TODO: Add optional waveform from NWB file
% if ...
%   spike.waveform = waveform;
% end
if isfield(spike, 'waveform')
    spike.dimord = '{chan}_lead_time_spike';
else
    spike.dimord = 'lead_time_spike'; % TODO: Or is it always {chan}_lead_time_spike?
end
if exist('nwb_trial_stime', 'var')
    spike.trial = trial;
    spike.trialtime = [nwb_trial_stime, nwb_trial_etime];
end
if exist('nwb_unit_info', 'var'), spike.nwb_unit_info = nwb_unit_info; end
if exist('nwb_stim_info', 'var'), spike.nwb_stim_info = nwb_stim_info; end
if exist('nwb_trial_info', 'var'), spike.nwb_trial_info = nwb_trial_info; end
if exist('nwb_elec_info', 'var'), spike.nwb_elec_info = nwb_elec_info; end



% Nested functions
function strct = constrained_set_to_struct(cset, desc)
    % Convert constrained set (often used in NWB) to Matlab structure (used
    % in FieldTrip). This comes in handy for e.g. stimulus or trial
    % information. 
    %
    % INPUT:    cset  = constrained set
    %           desc  = (optional) string denoting information about the
    %                   input set (e.g. for stimulus information from the NWB file
    %                   you might consider putting 'from nwb.stimulus_presentation')
    % OUTPUT:   strct = matlab structure containing information in set with
    %                   set keys now being fields in the structure. 
    
    if ~exist('desc', 'var')
        desc = '';
    end
    
    % Add stimulus presentation info field
    ks = cset.keys();

    % if key name contains a '.', convert to '_'
    k_field = cell(size(ks));
    for ik = 1:length(ks)
        k_field{ik} = strrep(ks{ik}, '.', '_');
    end

    % convert NWB set to matlab structure and add to ft structure
    strct = cell2struct(cell(length(k_field),1),k_field);
    for ik = 1:length(ks)
        k = ks(ik);
        if isa(cset.get(k).data, 'types.untyped.DataStub')
            strct.(k_field{ik}) = cset.get(k).data.load;
        elseif isa(cset.get(k).data, 'types.untyped.ObjectView')
            strct.(k_field{ik}) = {cset.get(k).data.path};
        end
    end
    strct.desc = desc;
end % eofun constrained_set_to_struct


end % eofun read_nwb



