function test_bug3184

% WALLTIME 00:10:00
% MEM 300mb

% TEST getdimord ft_datatype_source ft_selectdata

% it seems it already goes wrong at an earlier stage, since the data lacks
% a csdlabel field. This seems the consequence of an intermediate call to
% ft_selectdata, subselecting trials, which loses the csdlabel along the
% way. First look into this.
load(dccnpath(fullfile('/home/common/matlab/fieldtrip/data/test/latest/source/meg' , 'source_grid_mtmfft_fourier_trl_PCC_keepall_ctf151')));

assert(isfield(source.avg, 'csdlabel'));

cfg = [];
cfg.trials = 1:5;
source2 = ft_selectdata(cfg, source);
assert(isfield(source2, 'avg') && isfield(source2.avg, 'csdlabel') || isfield(source2, 'csdlabel'));

source3 = ft_datatype_source(source);
assert(isequal(getdimord(source3, 'mom'),'{pos}_ori_rpt'));
assert(isequal(getdimord(source3, 'csd'),'{pos}_ori_ori'));
assert(isequal(getdimord(source3, 'filter'),'{pos}_ori_chan'));

source.cumtapcnt = ones(5,1)*2;
source4 = ft_datatype_source(source);
assert(isequal(getdimord(source4, 'mom'),'{pos}_ori_rpttap'));
assert(isequal(getdimord(source4, 'csd'),'{pos}_ori_ori'));
assert(isequal(getdimord(source4, 'filter'),'{pos}_ori_chan'));

% also try a different data structure, obtained with mtmconvol
load(dccnpath(fullfile('/home/common/matlab/fieldtrip/data/test/latest/source/meg' , 'source_grid_mtmconvol_fourier_trl_PCC_keepall_ctf151')));

assert(isfield(source.avg, 'csdlabel'));

cfg = [];
cfg.trials = 1:5;
source2 = ft_selectdata(cfg, source);
assert(isfield(source2, 'avg') && isfield(source2.avg, 'csdlabel') || isfield(source2, 'csdlabel'));

source3 = ft_datatype_source(source);
assert(isequal(getdimord(source3, 'mom'),'{pos}_ori_rpt'));
assert(isequal(getdimord(source3, 'csd'),'{pos}_ori_ori'));
assert(isequal(getdimord(source3, 'filter'),'{pos}_ori_chan'));

source.cumtapcnt = ones(5,1)*2;
source4 = ft_datatype_source(source);
assert(isequal(getdimord(source4, 'mom'),'{pos}_ori_rpttap'));
assert(isequal(getdimord(source4, 'csd'),'{pos}_ori_ori'));
assert(isequal(getdimord(source4, 'filter'),'{pos}_ori_chan'));

% load the data that has been submitted with the bug report.
load(dccnpath(fullfile('/home/common/matlab/fieldtrip/data/test' , 'bug3184')));

%getdimord(ft_datatype_source(source_wvlt),'mom')

source_wvlt.avg.csdlabel = {'scandip'};
assert(isequal(getdimord(ft_datatype_source(source_wvlt),'mom'), '{pos}_ori_rpt'));




