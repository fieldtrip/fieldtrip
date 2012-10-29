function test_ft_apply_montage

% TEST test_ft_apply_montage
% TEST ft_apply_montage

pwdir = pwd;

cd('/home/common/matlab/fieldtrip/data/');
cfg = [];
cfg.dataset = 'Subject01.ds';
cfg.trl     = [1 1200 0];
cfg.continuous = 'yes';
data = ft_preprocessing(cfg);

mont          = [];
mont.tra      = -eye(151); %flip sign
mont.labelorg = ft_channelselection({'MEG'}, data.label);
mont.labelnew = ft_channelselection({'MEG'}, data.label);

data2 = ft_apply_montage(data, mont, 'keepunused', 'yes');
data3 = ft_apply_montage(data, mont, 'keepunused', 'no');
grad  = ft_apply_montage(data.grad, mont, 'keepunused', 'yes');

[a,b] = match_str(data2.label, data.label);
assert(all(all(data2.trial{1}(a(1:151),:)==-data.trial{1}(b(1:151),:)))==1);
assert(all(all(data.grad.tra(1:151,:)==-grad.tra(1:151,:)))==1);

cd(pwdir);
