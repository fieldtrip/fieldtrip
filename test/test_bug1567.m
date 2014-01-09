function test_bug1567

% MEM 1500mb
% WALLTIME 00:03:05

% TEST test_bug1567
fileloc = dccnfilename('/home/common/matlab/fieldtrip/data/test/bug1567');

datasets = {'test0001.eeg';
            'test0001.vhdr';
            'test0001.vmrk';
            'test0001-dat.dat';
            'test0001-dat.vhdr';
            'test0001-dat.vmrk';
            'test0002-seg.seg';
            'test0002-seg.vhdr';
            'test0002-seg.vmrk'};
%  load datasets

cfg = [];
cfg.continuous = 'yes';
for i=[1 2 4 5 7 8] % 1:numel(datasets) - exclude .vmrk 
   cfg.dataset = [fileloc filesep datasets{i}];
   if strcmp(datasets{i}(10:12),'seg')
     cfg.continuous = 'no';
   end;
   data = ft_preprocessing(cfg);
end;