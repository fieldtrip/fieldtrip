% test the ft_preprocessing for re-referencing to rest
% Li Dong $ 2020.3.3, (lidong@uestc.edu.cn)
clear all;
clc;

ft_defaults; 

load F:\GitHub_work\test_data\lf.mat % load leadfield calculated by fieldtrip
% load F:\GitHub_work\test_data\leadfieldMatrix.mat % load leadfield calculated by fieldtrip

cfg = [];
cfg.dataset     = 'F:\GitHub_work\test_data\test_sub01_avg.set';
cfg.reref       = 'yes';
cfg.refmethod     = 'rest';     %  if select 'rest','leadfield' is required.
cfg.leadfield = lf;            
%              The leadfield can be a matrix (channels X sources)
%              which is calculated by using the forward theory, based on
%              the electrode montage, head model and equivalent source
%              model. It can also be the output of ft_prepare_leadfield.m
%              (e.g. lf.leadfield) based on real head modal using FieldTrip.
% refchann = [];
% for i = 1:60
%     refchann{i} = num2str(i);
% end

cfg.refchannel = lf.label(1:60,1); % use first 60 channels
% cfg.refchannel     = {'all'};   % vector with indices of the selected channels 
%                               %  (re-referenced channels), or 'all'.
data_eeg        = ft_preprocessing(cfg);

