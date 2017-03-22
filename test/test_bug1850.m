function test_bug1850

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_prepare_neighbours ft_channelrepair
%
% http://bugzilla.fcdonders.nl/show_bug.cgi?id=1850

load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf275.mat'));

cfg=[];
cfg.method='template';
cfg.template='CTF275_neighb.mat';
n=ft_prepare_neighbours(cfg);

% get the 'full' list of channel names
for i=1:length(n)
    allchans{i,:}=n(i).label;
end

missingchans=setdiff(allchans,data.label);

% repair
cfg=[];
cfg.missingchannel=missingchans;
cfg.neighbours=n;
cfg.method='spline';
data_r=ft_channelrepair(cfg,data);
