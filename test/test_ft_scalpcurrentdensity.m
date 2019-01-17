function test_ft_scalpcurrentdensity

% MEM 7500mb
% WALLTIME 00:10:00

% TEST ft_sourceinterpolate ft_sourceplot

clear all
close all
addpath /project/3017042.02/fieldtrip;
ft_defaults;
%% All methods, with and without bad channels

% Generate some data
% data = [];
% data.trial{1} = randn(3,100);
% data.trial{2} = randn(3,100);
% data.time{1}  = (0:99)./100;
% data.time{2}  = (0:99)./100;
% data.label    = {'chan1';'chan2';'chan3'};

% cfg.elec.label = data_eeg.label;
% cfg.elec.elecpos = [1 1 1; 2 2 2; 3 3 3];
% cfg.elec.chanpos = [1 1 1; 2 2 2; 3 3 3];

% Load data 
load('/project/3017042.02/Log/EEG/postICA/MRIpav_001_postIC.mat');
load('/project/3017042.02/Log/EEG/EEG_electrode_pos/posLabels.mat');
elec = ft_read_sens('/project/3017042.02/Log/EEG/EEG_electrode_pos/MRIpav_001_001.pos');

% Duplicate data
baddata = data;
% Manipulate channel 1 to NaNs:
baddata.trial{1}(1,1) = NaN; % trial 1, channel 1, sample 1

% Loop over methods:
methods  = {'finite', 'spline'}; %, 'hjorth'}; % hjorth method needs neighbours
for i = 1:length(methods)
    cfg = [];
    cfg.method = string(methods(i));
    cfg.elec = elec;
    cfg.elec.label(1:67) =  posLabels; % overwrite to have channel labels (e.g. FPz) instead of numbers
    % Select channels also contained in data.label:
%     [dataindx, elecindx] = match_str(data.label, cfg.elec.label);
%     cfg.elec.chanpos = cfg.elec.chanpos(elecindx,:);
%     cfg.elec.chantype = cfg.elec.chantype(elecindx);
%     cfg.elec.chanunit = cfg.elec.chanunit(elecindx);
%     cfg.elec.elecpos = cfg.elec.elecpos(elecindx,:);
%     cfg.elec.label = cfg.elec.label(elecindx);  
    % Select neigbours for hjorth method:
%     tmpcfg = [];
%     tmpcfg.method = 'distance';
%     tmpcfg.neighbourdist = 1; % ???
%     tmpcfg.layout = 1; % ???
%     tmpcfg.channels = 'chan1'; % ???
%     % OR:
%     tmpcfg.elec = cfg.elec;
%     cfg.neighbours = ft_prepare_neighbours(tmpcfg, data);

    % a) Data without bad channels:
    fprintf('>>> START: Try out method %s without bad channels\n',string(methods(i)))
    scd = ft_scalpcurrentdensity(cfg,data);
    % b) Data with bad channels:
    fprintf('>>> START: Try out method %s with bad channels\n',string(methods(i)))
    scd = ft_scalpcurrentdensity(cfg,baddata);
    % any assertions to check that scd isn't just rubbish?
    if(sum(cellfun(@(C) any(isnan(C(:))), scd.trial))==0)
        fprintf('>>> No NaNs in output found, so all good\n')
    else
        fprintf('>>> NaNs found in output \n')        
    end
end

%END