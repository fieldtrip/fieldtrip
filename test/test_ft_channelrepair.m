function test_ft_channelrepair

% TEST test_ft_channelrepair
% TEST ft_channelrepair ft_datatype_sens fixsens ft_prepare_neighbours

datainfo = ref_datasets;

% get an MEG and an EEG set (hard-coded
eeginfo = datainfo(3);
meginfo = datainfo(7);

% do the MEG processing
fname = fullfile(meginfo.origdir,'latest', 'raw',meginfo.type,['preproc_' meginfo.datatype]);
load(fname);

cfg = [];
cfg.method = 'triangulation';
neighbours = ft_neighbourselection(cfg, data);

cfg = [];
cfg.badchannel = data.label(100);
cfg.neighbours = neighbours;
newdata = ft_channelrepair(cfg, data);

% % do the EEG processing: this does not work, because there's no example EEG data
% with sensor positions 
% fname = [eeginfo.origdir,'raw/',eeginfo.type,'preproc_',eeginfo.datatype];
% load(fname);
% 
% cfg = [];
% cfg.method = 'template';
% cfg.template = 'EEG1010_neighb.mat';
% neighbours = ft_neighbourselection(cfg);
% 
% cfg = [];
% cfg.badchannel = data.label(10);
% cfg.neighbours = neighbours;
% newdata = ft_channelrepair(cfg, data);


%% part 2 - missing channels and EEG data
% make use of bug941 data
% load data
if ispc
    home_dir = 'H:';
else    
    home_dir = '/home';
end
main_dir = fullfile(home_dir, 'common', 'matlab', 'fieldtrip', 'data', 'test');
bug_data = 'bug941.mat';
load(fullfile(main_dir, bug_data));

% treat as a bad channel
data_eeg_clean.elec = elec_new;
cfg = [];
cfg.badchannel = {'25'};
cfg.neighbours = neighbours;
data_eeg_repaired = ft_channelrepair(cfg,data_eeg_clean);

cfg = [];
cfg.channel = {'19','20','24','25','26'};
ft_databrowser(cfg, data_eeg_repaired);

% treat as a missing channel
data_eeg_miss = data_eeg_clean;
data_eeg_miss.label(25) = []; % remove channel 25

for i=1:numel(data_eeg_miss.trial) 
  data_eeg_miss.trial{i}(25, :) = []; % remove data from channel 25
end

cfg = [];
cfg.missingchannel = {'25'};
cfg.neighbours = neighbours;
data_eeg_interp = ft_channelrepair(cfg,data_eeg_miss);

cfg = [];
cfg.channel = {'19','20','24','25','26'};
ft_databrowser(cfg, data_eeg_interp);

% check for each trial whether the value for channel 25 is in between its 
% neighbours (it's now the last channel), and equals data_eeg_repaired
cfg.neighbours = [19 20 24 26];
for tr=1:numel(data_eeg_interp.trial)
  tmp = repmat(data_eeg_interp.trial{tr}(end, :), 4, 1);
  if all(tmp < data_eeg_interp.trial{tr}(cfg.neighbours, :)) | ...
      all(tmp > data_eeg_interp.trial{tr}(cfg.neighbours, :))
    error(['The average is not in between its channel neighbours at for trial ' num2str(tr)]);
  elseif ~all(data_eeg_interp.trial{tr}(end, :) == data_eeg_repaired.trial{tr}(25, :))
    error('The reconstruction of the same channel differs when being treated as a missing channel compared to a bad channel');
  else
    fprintf('trial %i is fine\n', tr);
  end
end

%% part 3 - spherical spline interpolation
% again make use of bug941 data
% load data

% treat as a bad channel
cfg = [];
cfg.badchannel = {'25'};
cfg.neighbours = neighbours;
cfg.method     = 'spline';
data_eeg_repaired_spline = ft_channelrepair(cfg,data_eeg_clean);

cfg = [];
cfg.channel = {'19','20','24','25','26'};
ft_databrowser(cfg, data_eeg_repaired_spline);

% treat as a missing channel
data_eeg_miss = data_eeg_clean;
data_eeg_miss.label(25) = []; % remove channel 25

for i=1:numel(data_eeg_miss.trial) 
  data_eeg_miss.trial{i}(25, :) = []; % remove data from channel 25
end

cfg = [];
cfg.missingchannel = {'25'};
cfg.neighbours = neighbours;
cfg.method     = 'spline';
data_eeg_interp_spline = ft_channelrepair(cfg,data_eeg_miss);

cfg = [];
cfg.channel = {'19','20','24','25','26'};
ft_databrowser(cfg, data_eeg_interp_spline);

% check for each trial whether the value for channel 25 is in between its 
% neighbours (it's now the last channel), and equals data_eeg_repaired
cfg.neighbours = [19 20 24 26];
for tr=1:numel(data_eeg_interp_spline.trial)
  %
  % spline interpolation does not yield to an average of all channels!
  %
  
  %tmp = repmat(data_eeg_interp.trial{tr}(end, :), 4, 1);  
  %if all(tmp < data_eeg_interp.trial{tr}(cfg.neighbours, :)) | ...
  %    all(tmp > data_eeg_interp.trial{tr}(cfg.neighbours, :))
  %  error(['The average is not in between its channel neighbours at for trial ' num2str(tr)]);
  %else
   if ~mean(abs(data_eeg_interp_spline.trial{tr}(end, :) - data_eeg_repaired_spline.trial{tr}(25, :))) < 10e6*eps % i.e. nearly ==0
    disp(['mean error is: ' num2str(mean(abs(data_eeg_interp_spline.trial{tr}(end, :) - data_eeg_repaired_spline.trial{tr}(25, :))))]); 
    error('The reconstruction of the same channel differs when being treated as a missing channel compared to a bad channel');
  else
    fprintf('trial %i is fine\n', tr);
  end
end