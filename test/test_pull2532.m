function test_pull2532

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_componentanalysis ft_rejectcomponent ft_apply_montage
% DATA none

% test the rank of the data after PCA and backprojection
% test the size of elecpos and chanpos after PCA and backprojection

%%

nchan = 10;
nsample = 1000;
fsample = 1000;

data = [];
data.label = {};
for i=1:nchan
  data.label{i} = sprintf('eeg%02d', i);
end
data.time{1} = (1:nsample)/fsample;
data.trial{1} = randn(nchan, nsample);
for i=1:nchan
  % channel 1, is the largest, followed by channel 2, etc.
  data.trial{1}(i,:) = data.trial{1}(i,:)/i;
end

elec = [];
elec.label = data.label;
elec.elecpos = zeros(nchan,3);
elec.unit = 'cm';
for i=1:nchan
  elec.elecpos(i,3) = i;
end
elec.chanpos = elec.elecpos;

% add the electrode structure
% weird stuff will happen with the chanpos after the comp+invcomp projection
data.elec = elec;

%%

cfg = [];
cfg.method = 'pca';
comp = ft_componentanalysis(cfg, data);
assert(numel(comp.label)==nchan);
assert(rank(comp.trial{1})==nchan);
assert(size(comp.elec.elecpos,1)==nchan);
assert(size(comp.elec.chanpos,1)==nchan);

cfg = [];
cfg.component = [];
clean = ft_rejectcomponent(cfg, comp);
assert(numel(clean.label)==nchan);
assert(rank(clean.trial{1})==nchan);
assert(size(comp.elec.elecpos,1)==nchan);
assert(size(comp.elec.chanpos,1)==nchan);

cfg = [];
cfg.component = 6:10; % remove 5 components
clean = ft_rejectcomponent(cfg, comp);
assert(numel(clean.label)==nchan);
assert(rank(clean.trial{1})==5);
assert(size(comp.elec.elecpos,1)==nchan);
assert(size(comp.elec.chanpos,1)==nchan);

%%

cfg = [];
cfg.method = 'pca';
cfg.numcomponent = 5;
comp = ft_componentanalysis(cfg, data);
assert(numel(comp.label)==5);
assert(rank(comp.trial{1})==5);
assert(size(comp.elec.elecpos,1)==nchan);
assert(size(comp.elec.chanpos,1)==5); % reduced

cfg = [];
cfg.component = [];
clean = ft_rejectcomponent(cfg, comp);
assert(numel(clean.label)==nchan);
assert(rank(clean.trial{1})==5); % there should be 5 left
assert(size(clean.elec.elecpos,1)==nchan);
assert(size(clean.elec.chanpos,1)==nchan); % back to nchan

cfg = [];
cfg.component = [3 4 5]; % remove another three
clean = ft_rejectcomponent(cfg, comp);
assert(numel(clean.label)==nchan);
assert(rank(clean.trial{1})==2); % there should be 2 left
assert(size(clean.elec.elecpos,1)==nchan);
assert(size(clean.elec.chanpos,1)==nchan); % back to nchan

%%

cfg = [];
cfg.method = 'pca';
cfg.numcomponent = 5;
comp = ft_componentanalysis(cfg, data);
assert(numel(comp.label)==5);
assert(rank(comp.trial{1})==5);
assert(size(comp.elec.elecpos,1)==nchan);
assert(size(comp.elec.chanpos,1)==5); % reduced

cfg = [];
cfg.component = 1; % remove 1 out of 5
clean = ft_rejectcomponent(cfg, comp, data); % now on the original data
assert(numel(clean.label)==nchan);
assert(rank(clean.trial{1})==9); % there should be 9 left
assert(size(clean.elec.elecpos,1)==nchan);
assert(size(clean.elec.chanpos,1)==nchan); % back to nchan
