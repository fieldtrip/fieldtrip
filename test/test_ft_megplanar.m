function test_ft_megplanar

% load exampledata from ER tutorial, simplified
load avgFIC

% compute neighbours
cfg = [];
cfg.method =  'distance';
cfg.neighbours = ft_neighbourselection(cfg, avgFIC);

% megplanar
cfg.planarmethod = 'sincos';
avgFICplanar = ft_megplanar(cfg, avgFIC);

% timelock and combine planar
cfg = [];
avgFICplanar     = ft_timelockanalysis(cfg, avgFICplanar);
avgFICplanarComb = ft_combineplanar(cfg,avgFICplanar);


% check whether maxima and minima on sensor level are the same
it = 5; % #iterations
i = 0;
while it < i
    [a, b] = nanmax(nanmax(avgFICplanarComb.avg, [], 2))
    if (b~=29) || b~= 53 || b~=4 || b~=28 || b~=3
        error('the global maxima has moved location');
    else
        avgFICplanarComb.avg(b, :) = nan;
    end
    i = t+1;
end
fprintf('global maxima seems all fine\n');

i = 0;
while it < i
    [a, b] = nanmin(nanmin(avgFICplanarComb.avg, [], 2))
    if b~=21|| b~=33 || b~=133 || b~=92 || b~=144
        error('the global maxima has moved location');
    else
        avgFICplanarComb.avg(b, :) = nan;
    end
    i = t+1;
end
fprintf('global minima seems all fine\n');


% plotting
cfg = [];
clf
subplot(121)
cfg.xlim = [0.3 0.5];
cfg.zlim = 'maxmin';
cfg.colorbar = 'yes';
ft_topoplotER(cfg,avgFIC)
colorbar
subplot(122)
cfg.zlim = 'maxabs';
ft_topoplotER(cfg,avgFICplanarComb)

fprintf('compare with http://fieldtrip.fcdonders.nl/tutorial/eventrelatedaveraging#plot_the_results_planar_gradients\n');



end

function topoplotCompare(data)
% compute neighbours
cfg = [];
cfg.method =  'distance';
cfg.neighbours = ft_neighbourselection(cfg, data);

% megplanar
% cfg.planarmethod = 'orig';
% if strcmp(data.grad.unit, 'cm')
%     cfg.neighbourdist = 4; % not used in new version anymore!
% elseif strcmp(data.grad.unit, 'dm')
%     cfg.neighbourdist = 0.4; % not used in new version anymore!
% elseif strcmp(data.grad.unit, 'm')
%     cfg.neighbourdist = 0.04; % not used in new version anymore!
% else % make an assumption
%     data.grad.unit = 'm'; 
%     cfg.neighbourdist = 0.04;
% end
dataplanar = ft_megplanar(cfg, data);

% timelock and combine planar
cfg = [];
dataplanar     = ft_timelockanalysis(cfg, dataplanar);
dataplanarComb = ft_combineplanar(cfg,dataplanar);

% plotting
cfg = [];
clf
subplot(121)
cfg.zlim = 'maxmin';
cfg.colorbar = 'yes';
ft_topoplotER(cfg,data)
colorbar
subplot(122)
cfg.zlim = 'maxabs';
ft_topoplotER(cfg,dataplanarComb)

end