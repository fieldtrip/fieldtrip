function test_issue2585

filename = dccnpath('/project/3031000.02/test/test_issue2585.mat');

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY fixbalance
% DATA private


% FieldTrip test for projectors, original contribution by rcassani, to
% demonstrate the issue

% %% Preprocessing 1
% cfg            = [];
% cfg.dataset    = 'ArtifactMEG.ds';
% cfg.continuous = 'yes'; % see https://www.fieldtriptoolbox.org/faq/continuous/
% cfg.channel    = {'MEG', 'MEGREF', 'EOG', 'ECG', 'jaw', 'neck'};
% data = ft_preprocessing(cfg);
% % Resample
% cfg              = [];
% cfg.resamplefs   = 300;
% cfg.detrend      = 'no';
% data = ft_resampledata(cfg, data);
% 
% 
% %% ICA decomposition
% 
% % ICA decomposition 1
% cfg              = [];
% cfg.method       = 'pca';
% cfg.numcomponent = 20;
% cfg.channel      = {'MEG', 'MEGREF'};
% data_comp_1 = ft_componentanalysis(cfg, data); % using the data without atypical artifacts
% % Remove ICs
% cfg = [];
% cfg.component = [10 12]; % to be removed
% data_clean_1 = ft_rejectcomponent(cfg, data_comp_1);
% 
% % ICA decomposition 2
% cfg              = [];
% cfg.method       = 'pca';
% cfg.numcomponent = 15;
% cfg.channel      = {'MEG', 'MEGREF'};
% data_comp_2 = ft_componentanalysis(cfg, data_clean_1); % using the data without atypical artifacts
% % Remove IC
% cfg = [];
% cfg.component = [9]; % to be removed
% data_clean_2 = ft_rejectcomponent(cfg, data_comp_2);
% 
% %% MONTAGE ORDER
% grad = data_clean_2.grad;
% if isfield(grad.balance, 'previous')
%     % Convert to new sturecture using code from `fixbalance.m`
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if ischar(grad.balance.current)
%       grad.balance.current = {grad.balance.current};
%     end
% 
%     if isfield(grad.balance, 'previous')
%       % concatenate them and keep a single list
%       grad.balance.current = cat(1, grad.balance.previous(:), grad.balance.current(:))';
%       grad.balance = rmfield(grad.balance, 'previous');
%     end
% 
%     % remove the 'none'
%     grad.balance.current = grad.balance.current(~strcmp(grad.balance.current, 'none'));
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end
% 
% grad.balance.current
% 
% % After updates in September 2025
% % grad.balance.current{1} is the first and 
% % grad.balance.current{end} is the last applied montage
% % grad.balance.current = {'comp'}    {'invcomp'}    {'comp1'}    {'invcomp1'}
% 
% % Before updates in September 2025
% % grad.balance.current       was the last applied montage,
% % grad.balance.previous{1}   was the second-last applied montage and
% % grad.balance.previous{end} was the first applied montage
% % grad.balance.current = {'comp2'}    {'invcomp1'}    {'comp1'}    {'invcomp2'}
% %
% % `fixbalance.m` needs to reverse grad.balance.previous  before concatenate

% the file test_issue2585 contains 2 grad structures:
%  - gradnew, as obtained directly from a sequence of calls to componentanalysis/rejectcomponent with the new grad representation
%  - grad, as obtained from an older version of fieldtrip
%
% NOTE: THERE ARE NUMERIC DIFFERENCES IN THE RESULTING TRA MATRIX, WHICH I
% (JM) SUSPECT TO BE THE CONSEQUENCE OF SOMEHOW SLIGHT UNCONTROLLED
% DIFFERENCES IN THE CODEBASE, DIFFERENCE IN CHANNEL/COIL ORDER?

[ver, ftpath] = ft_version;
cd(fullfile(ftpath, 'private'));

load(filename);
gradfixed = fixbalance(grad);
assert(isequal(gradfixed.balance.current, {'comp1' 'invcomp1' 'comp2' 'invcomp2'}));
assert(isalmostequal(gradnew.tra,gradfixed.tra,'reltol',0.00001))
