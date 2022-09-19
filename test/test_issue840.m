function test_issue840

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_megplanar ft_combineplanar ft_appenddata

%%

if false
  % let's make three small data structures
  % the data is from https://openneuro.org/datasets/ds000117/versions/1.0.3
  
  filename = {
    'sub-01_ses-meg_task-facerecognition_run-01_meg.fif'
    'sub-01_ses-meg_task-facerecognition_run-02_meg.fif'
    'sub-01_ses-meg_task-facerecognition_run-03_meg.fif'
    };
  
  
  for i=1:numel(filename)
    cfg = [];
    cfg.trialfun = 'ft_trialfun_general';
    cfg.trialdef.triallength = 1;
    cfg.trialdef.ntrials = 10;
    cfg.dataset = filename{i};
    cfg = ft_definetrial(cfg);
    data{i} = ft_preprocessing(cfg);
  end
  
  data1 = data{1};
  data2 = data{2};
  data3 = data{3};
  clear data
  save data1 data1
  save data2 data2
  save data3 data3
  
else
  % load it from the mat files
  load(dccnpath('/home/common/matlab/fieldtrip/data/test/issue840/data1.mat'));
  load(dccnpath('/home/common/matlab/fieldtrip/data/test/issue840/data2.mat'));
  load(dccnpath('/home/common/matlab/fieldtrip/data/test/issue840/data3.mat'));
  
end

%%

cfg = [];
data1_combined = ft_combineplanar(cfg, data1);
data2_combined = ft_combineplanar(cfg, data2);
data3_combined = ft_combineplanar(cfg, data3);

cfg = [];
data_combined_append = ft_appenddata(cfg, data1_combined, data2_combined, data3_combined);
assert( isfield(data_combined_append, 'elec')); % this one is consistent
assert(~isfield(data_combined_append, 'grad')); % this one is not consistent

%%

cfg = [];
data_append = ft_appenddata(cfg, data1, data2, data3);
assert( isfield(data_append, 'elec')); % this one is consistent
assert(~isfield(data_append, 'grad')); % this one is not consistent

%%

cfg = [];
cfg.channel = 'MEG*3';
cfg.method = 'distance';
cfg.grad = ft_convert_units(data1.grad, 'mm');
cfg.neighbourdist = 40;
neighbours = ft_prepare_neighbours(cfg);

%%

cfg = [];
cfg.channel = 'MEG*3';
cfg.neighbours = neighbours;

% the FT_MEGPLANAR function refuses to work on neuromag306 data. Here we carefully
% selected only the magnetometers (both as cfg.channel and in cfg.neighbours) so in
% principle it should work fine. The hack I am using is to make FT_MEGPLANAR think it
% is an itab153 sensor array
data1.grad.type = 'itab153';
data2.grad.type = 'itab153';
data3.grad.type = 'itab153';

data1_planar = ft_megplanar(cfg, data1);
data2_planar = ft_megplanar(cfg, data2);
data3_planar = ft_megplanar(cfg, data3);

%%

cfg = [];
data_planar_append = ft_appenddata(cfg, data1_planar, data2_planar, data3_planar);
assert( isfield(data_planar_append, 'elec')); % this one is consistent
assert(~isfield(data_planar_append, 'grad')); % this one is not consistent

%%

cfg = [];
cfg.channel = {'*_dV', '*_dV'};
data1_planar_combined = ft_combineplanar(cfg, data1_planar);
data2_planar_combined = ft_combineplanar(cfg, data2_planar);
data3_planar_combined = ft_combineplanar(cfg, data3_planar);

%%

cfg = [];
data_planar_combined_append = ft_appenddata(cfg, data1_planar_combined, data2_planar_combined, data3_planar_combined);
assert( isfield(data_planar_combined_append, 'elec')); % this one is consistent
assert(~isfield(data_planar_combined_append, 'grad')); % this one is not consistent

%%

% this fails because it cannot determine the sensor type: the data contains elec, but no grad

% However, if it were to contain grad, it would probably also fail since I modified
% the grad.type above from neuromag306 into itab153. Overall, this is a pretty crappy
% test script. I'll keep it nevertheless, since there are some usefiel parts in it.

if false
  cfg = [];
  cfg.channel = {'*_dV', '*_dV'};
  data_planar_append_combined = ft_combineplanar(cfg, data_planar_append);
end

% Error using ft_checkdata (line 573)
% This function requires data with an 'ctf151_planar', 'ctf275_planar', 'neuromag122',
% 'neuromag306', 'bti248_planar', 'bti148_planar', 'itab153_planar', 'yokogawa160_planar',
% 'yokogawa64_planar' or 'yokogawa440_planar' sensor array.




