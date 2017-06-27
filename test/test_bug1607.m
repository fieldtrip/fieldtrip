function test_bug1607

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_channelrepair ft_topoplotER

% During preprocessing I lost some channels which I got back through
% ft_channelrepair after which the channel order changed per subject.
% I realigned the channel order and the associated data structure of each subject
% to one reference order after which the topoplots look completely different.
% They look even worse than before realignment. Anybody got a clue where the
% confusion arises?
% Attached a script and two data files, one with reference data and one with the
% to be aligned data.

% solution
% was not a bug, instead jonas was indexing the wrong way around:
%
%     nERPdata_nD_left{isubject}.avg(:,:) = ERPdata_nD_left{isubject}.avg(loc,:);     % realign data to ref channel order
% instead of
%     nERPdata_nD_left{isubject}.avg(loc,:) = ERPdata_nD_left{isubject}.avg(:,:);     % realign data to ref channel order

ERPdata_nD_left = [];
reference_labels = [];

%load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1607/06_control_ICA_clean.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1607/ERPdata.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1607/reference.mat'));

% reference_labels = ICA_clean.label;

nERPdata_nD_left = ERPdata_nD_left;
% nERPdata_D_left = ERPdata_D_left;
% nERPdata_nD_right = ERPdata_nD_right;
% nERPdata_D_right = ERPdata_D_right;

for isubject = 1:1%length(ERPdata_nD_left)
%     disp(num2str(isubject))
    [a loc] = ismember(ERPdata_nD_left{isubject}.label, reference_labels);  % compare the channel labeling of each subject with a reference
%     for isens = 1 : 64
%         if loc(isens) ~= isens
%             disp(['changing subj ' num2str(isubject) ' sensnr ' num2str(isens)]);
%             break;
%         end
%     end
%     disp(num2str(any(find((ERPdata_nD_left{isubject}.avg == ERPdata_nD_left{isubject}.avg(loc,:))==0))));    
    nERPdata_nD_left{isubject}.label = reference_labels;    % relabel the channels according to the reference
%     nERPdata_nD_left{isubject}.avg(:,:) = ERPdata_nD_left{isubject}.avg(loc,:);     % realign data to ref channel order
    nERPdata_nD_left{isubject}.avg(loc,:) = ERPdata_nD_left{isubject}.avg(:,:);     % realign data to ref channel order
    nERPdata_nD_left{isubject}.dof(:,:) = ERPdata_nD_left{isubject}.dof(loc,:);     % realign data to ref channel order

end
% after realigning the data the plots look completely different, compare
% old and new data

cfg = [];
cfg.layout      = 'biosemi64.lay';
cfg.xlim        = [.04 .06];
cfg.channel     = 'all';
cfg.interactive = 'no';
cfg.highlight   = 'yes';
cfg.highlightchannel = reference_labels(end-6:end);

figure; ft_topoplotER(cfg,ERPdata_nD_left{1});      % plot original data 
figure; ft_topoplotER(cfg,nERPdata_nD_left{1});     % plot corrected data
