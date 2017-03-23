function test_old_appenddata_2datasets

% MEM 1gb
% WALLTIME 00:10:00


% this script tests ft_appenddata when the input is obtained from 2
% different datafiles

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/ctf275/A0132_Aud-Obj-Recognition_20051115_02.ds'));
headerfile = 'A0132_Aud-Obj-Recognition_20051115_02.res4';
datafile   = 'A0132_Aud-Obj-Recognition_20051115_02.meg4';
%hdr        = ft_read_header(headerfile);

cfg          = [];
cfg.datafile = datafile;
cfg.trl      = [[1001:1000:10001]' [2000:1000:11000]' round(randn(10,1)*100)];
cfg.trl(:,4) = [ones(5,1); ones(5,1)*2];
cfg.continuous = 'yes';
data1          = ft_preprocessing(cfg);
cfg.trl(:,1:2) = cfg.trl(:,1:2) + 11000;
cfg.trl(:,3)   = round(randn(10,1)*100);
data2          = ft_preprocessing(cfg);

data3a             = ft_appenddata([], data1, data2);
data1.cfg.datafile = 'file1';
data2.cfg.datafile = 'file2';
data3b             = ft_appenddata([], data1, data2);
