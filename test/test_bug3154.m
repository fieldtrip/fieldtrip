function test_bug3154

% WALLTIME 00:20:00
% MEM 4gb

%%

if true
  load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/sensor_analysis/subjectK.mat'));
  data = data_left;
  nMEG = 151;
  senstype = 'ctf151';
else
  load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/natmeg/preprocessing/data_MEG.mat'));
  data = data_MEG;
  nMEG = 306;
  senstype = 'neuromag306';
end

%%

type1 = ft_chantype(data);
type2 = ft_chantype(data.label);
type3 = ft_chantype(data.hdr);
type4 = ft_chantype(data.hdr.label);
type5 = ft_chantype(data.grad);
type6 = ft_chantype(data.grad.label);

assert(mean(strcmp(type1, 'unknown'))<0.5); % this one failed on 2016-06-27
assert(mean(strcmp(type2, 'unknown'))<0.5);
assert(mean(strcmp(type3, 'unknown'))<0.5);
assert(mean(strcmp(type4, 'unknown'))<0.5);
assert(mean(strcmp(type5, 'unknown'))<0.5);
assert(mean(strcmp(type6, 'unknown'))<0.5);

%%

sel2 = ft_channelselection('MEG', data.label);
sel4 = ft_channelselection('MEG', data.hdr.label);
sel6 = ft_channelselection('MEG', data.grad.label);

assert(numel(sel2)==nMEG);
assert(numel(sel4)==nMEG);
assert(numel(sel6)==nMEG);

%%
% the neuromag306 data was incorrectly detected as neuromag306_combined

assert(strcmp(ft_senstype(data),            senstype));
assert(strcmp(ft_senstype(data.label),      senstype));
assert(strcmp(ft_senstype(data.hdr),        senstype));
assert(strcmp(ft_senstype(data.hdr.label),  senstype));
assert(strcmp(ft_senstype(data.grad),       senstype));
assert(strcmp(ft_senstype(data.grad.label), senstype));

%%

cfg = [];
cfg.channel = 'MEG';

data_sel = ft_selectdata(cfg, data);
assert(numel(data_sel.label)==nMEG);

data_sel = ft_preprocessing(cfg, data);
assert(numel(data_sel.label)==nMEG);

