function test_bug2133

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug2133
% TEST ft_read_header read_eeg_mff

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeg/egi/NS500Sine6Hz.mff');

global testvariable

testvariable = 1;

hdr = ft_read_header(filename,'dataformat','egi_mff_v2', 'headerformat','egi_mff_v2');

assert(exist('testvariable', 'var')==1);

% The MATLAB isglobal function is obsolete and will be discontinued in a future version of MATLAB.
% Hence we need a local replication of this functionality
list = whos;
name = {list.name};
glob = [list.global];
x = any(strcmp('testvariable', name(glob)));
assert(x);

% it used to read
% assert(isglobal('testvariable'));

assert(testvariable==1);
