function test_bug941

% MEM 1500mb
% WALLTIME 00:03:08

% TEST test_bug941
% TEST ft_channelrepair ft_databrowser

% load data
if ispc
    home_dir = 'H:';
else    
    home_dir = 'home';cd /home/common/matlab/fieldtrip/data/test/
end
main_dir = fullfile(home_dir, 'common', 'matlab', 'fieldtrip', 'data', 'test');
bug_data = 'bug941.mat';
load(fullfile(main_dir, bug_data));

% code that produces the error
data_eeg_clean.elec = elec_new;
cfg = [];
cfg.badchannel = {'25'};
cfg.neighbours = neighbours;
data_eeg_repaired = ft_channelrepair(cfg,data_eeg_clean);

cfg = [];
cfg.channel = {'19','20','24','25','26'};
ft_databrowser(cfg, data_eeg_repaired);

% check for each trial whether the value for channel 25 is in between its 
% neighbours
cfg.neighbours = [19 20 24 26];
for tr=1:numel(data_eeg_repaired.trial)    
    ishigher = false;
    islower = false;
    tmp = repmat(data_eeg_repaired.trial{tr}(25, :), 4, 1);
    for t=1:size(data_eeg_repaired.trial{tr}, 2)        
        if (all(tmp(:, t) > data_eeg_repaired.trial{tr}([19 20 24 26], t)) || ...
            all(tmp(:, t) < data_eeg_repaired.trial{tr}([19 20 24 26], t)))            
            error(['The average is not in between its channel neighbours at time point ' num2str(data_eeg_repaired.time{tr}(t)) ' for trial ' num2str(tr)]);
        end
    end    
end

