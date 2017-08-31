function test_bug3337

% WALLTIME 00:10:00
% MEM 1gb

%%
data = [];
data.label = {'1', '2'};
for i=1:3
  data.trial{i} = randn(2, 1000);
  data.time{i} = (1:1000)/1000;
  data.sampleinfo(i,:) = [(i-1)*1000+1 i*1000];
end

data.trial{1}(:,30:40)  = nan;
data.trial{1}(:,80:90)  = nan; % 2nd segment within one trial
data.trial{2}(1,100)    = nan; % one sample
data.trial{3}(1,:)      = nan; % whole trial

cfg = [];
[cfg, artifact] = ft_artifact_nan(cfg, data);

assert(size(artifact,1)==4)
assert(size(artifact,2)==2)

correct = [
    30          40
    80          90
  1100        1100
  2001        3000
  ];

assert(isequal(artifact, correct));

%%

cfg.artfctdef.reject = 'nan';
datanan = ft_rejectartifact(cfg, data);
assert(isequaln(data.trial{1}, datanan.trial{1}));

%%

cfg.artfctdef.minaccepttim = 0;
cfg.artfctdef.reject = 'partial';
datapartial = ft_rejectartifact(cfg, data);
assert(numel(datapartial.trial)==5);
