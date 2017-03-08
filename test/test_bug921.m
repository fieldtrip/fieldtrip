function test_bug921

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_selectdata ft_selectdata_old ft_selectdata_new

% See also bug 798 that was reported by Yoni, from which I am reusing the data to test this bug

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug798'));
load t2_subj1.mat
load t2_subj1_null
load t2_subj2
load t2_subj2_null

argin{1} = ft_checkdata(t2_subj1, 'datatype', 'freq');
argin{2} = ft_checkdata(t2_subj1_null, 'datatype', 'freq');
argin{3} = ft_checkdata(t2_subj2, 'datatype', 'freq');
argin{4} = ft_checkdata(t2_subj2_null, 'datatype', 'freq');

% the following should select and concatenate only the powspctrm field
% but not the prob and stat fields
data = ft_selectdata(argin{:}, 'param', 'powspctrm');

assert(~isfield(data, 'stat'));
assert(~isfield(data, 'prob'));
assert(strcmp(data.dimord, 'rpt_chan_freq'));

assert(size(data.powspctrm,1)==4);
assert(size(data.powspctrm,2)==274);
assert(size(data.powspctrm,3)==1);
