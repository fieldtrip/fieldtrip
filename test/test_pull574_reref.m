function test_pull574_reref

% MEM 2gb
% WALLTIME 01:00:00
% DEPENDENCY ft_preprocessing preproc ft_preproc_reference

% test the ft_preprocessing for re-referencing to rest
% Li Dong $ 2020.3.10, (Lidong@uestc.edu.cn) and J.M.Schoffelen

load(dccnpath('/home/common/matlab/fieldtrip/data/test/pull574.mat')); % load data and leadfield calculated by FieldTrip

% -------------------------------------------------------------------------
% A snippet of code to calculate leadfield using FieldTrip
% cfg                        = [];
% cfg.elec                   = elec_aligned; % aligned eletrode positions
% cfg.grid                   = grid; % dipole positions
% cfg.headmodel              = vol;  % headmodel created by 'ft_prepare_headmodel'
% cfg.grid.resolution        = 6; % number for automatic sourcemodel generation
% cfg.normalize       = 'yes';
% [lf]        = ft_prepare_leadfield(cfg); % calculate leadfield
% -----------
% An full workflow of creating leadfield using real headmodel (BEM and FEM)
% can be seen in:
% http://www.fieldtriptoolbox.org/workshop/baci2017/forwardproblem/

% A full example script for leadfield based on FEM headmodel can be seen in
% http://www.fieldtriptoolbox.org/development/project/example_fem/
% -------------------------------------------------------------------------
% rest re-referencing
cfg           = [];
cfg.reref     = 'yes';
cfg.refmethod = 'rest';     %  if select 'rest','leadfield' is required.
cfg.leadfield = lf;
%              The leadfield can be a matrix (channels X sources)
%              which is calculated by using the forward theory, based on
%              the electrode montage, head model and equivalent source
%              model. It can also be the output of ft_prepare_leadfield.m
%              (e.g. lf.leadfield) based on real head modal using FieldTrip.

% cfg.refchannel = data.label([1:3,5:60],1); % use first 60 channels
cfg.refchannel = {'all'};    % vector with indices of the selected channels
%                                  %  (re-referenced channels), or 'all'.
data_eeg_rest  = ft_preprocessing(cfg,data);
% -------------------------------------------------------------------------
% median re-referencing
cfg = [];
cfg.reref       = 'yes';
cfg.refmethod   = 'median';    % median
cfg.refchannel  = {'all'};     % use first 60 channels

data_eeg_median = ft_preprocessing(cfg,data);

% plot the results

figure;
k = 1; % kth channel
data1 = data_eeg_rest.trial{1,1};
data2 = data_eeg_median.trial{1,1};

plot(data1(k,1:100),'-b');
hold on;
plot(data2(k,1:100),'-r');
plot(data1(k,1:100)-data2(k,1:100),'-.k')
hold off;grid on;
legend('rest','median','difference');
