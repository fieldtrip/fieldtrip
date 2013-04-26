function test_bug2133

% TEST test_bug2133
% TEST ft_read_header read_eeg_mff

% translate beween windows h:\ and linux /home path
filename = dccnfilename('/home/common/matlab/fieldtrip/data/test/original/eeg/egi/NS500Sine6Hz.mff');

global testvariable

testvariable = 1;

hdr = ft_read_header('SE7_T7_S5_GAV.mff','dataformat','egi_mff_v2', 'headerformat','egi_mff_v2');

assert(exist(testvariable, 'var'));
assert(isglobal(testvariable));
assert(testvariable==1);

