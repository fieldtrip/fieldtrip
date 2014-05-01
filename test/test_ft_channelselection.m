function test_ft_channelselection

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_ft_channelselection 
% TEST ft_senstype ref_datasets

datasets = ref_datasets;

for i=1:size(datasets,2);
    filename = fullfile(datasets(i).origdir, 'original', datasets(i).type, datasets(i).datatype, datasets(i).filename);
    
    if isempty(datasets(i).dataformat)
        hdr = ft_read_header(filename);
    else
        hdr = ft_read_header(filename,'headerformat',datasets(i).dataformat);
    end
    
    if ~ft_senstype(hdr, datasets(i).senstype)
        error(['incorrect senstype detection: ' datasets(i).datatype]);
    else
        type=ft_senstype(hdr);
    end
    
    if ~isnan(datasets(i).numeeg) && length(ft_channelselection('EEG', hdr.label, type))~=datasets(i).numeeg
        error('incorrect number of EEG channels');
    end
    if ~isnan(datasets(i).nummeg) && length(ft_channelselection('MEG', hdr.label, type))~=datasets(i).nummeg
        error('incorrect number of MEG channels');
    end
    if ~isnan(datasets(i).numeog) && length(ft_channelselection('EOG', hdr.label, type))~=datasets(i).numeog
        error('incorrect number of EOG channels');
    end
    if ~isnan(datasets(i).numecg) && length(ft_channelselection('ECG', hdr.label, type))~=datasets(i).numecg
        error('incorrect number of ECG channels');
    end
    if ~isnan(datasets(i).numemg) && length(ft_channelselection('EMG', hdr.label, type))~=datasets(i).numemg
        error('incorrect number of EMG channels');
    end

    %exceptions for 'neuromag306' and 'yokogawa160' MEG datasets:
    %selecting gradiometer/magnetometer sensors
    if strcmp(type,'neuromag306') || strcmp(type,'yokogawa160')
        if ~isnan(datasets(i).numemg) && length(ft_channelselection('MEGMAG', hdr.label, type))~=datasets(i).nummegmag
            error('incorrect number of MEGMAG channels');
        end
        if ~isnan(datasets(i).numemg) && length(ft_channelselection('MEGGRAD', hdr.label, type))~=datasets(i).nummeggrad
            error('incorrect number of MEGGRAD channels');
        end
    end
end
