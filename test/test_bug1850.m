function test_bug1850

% TEST test_bug1850
% TEST ft_prepare_neighbours ft_channelrepair
%
% http://bugzilla.fcdonders.nl/show_bug.cgi?id=1850


load test_bug1850; % load some CTF data
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