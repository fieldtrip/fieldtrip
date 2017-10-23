function test_bug2460

% WALLTIME 00:10:00
% MEM 1gb

% TEST ft_checkdata

% create pseudo dataset
data = [];
data.stat  = rand(4,5,6)*3;
data.prob  = rand(4,5,6);

% data.dim   = [4 5 6]; % this makes getdimord happy

% create two arbitrary mask fields
data.mask    = data.stat>1;
data.mask1   = data.stat>2;

data.inside  = 1:4*5*6/2;
data.outside = 4*5*6/2+1:4*5*6;
data.pos     = rand(4*5*6,3);
data.unit    = 'cm';

% this is how ft_checkdata is called from ft_sourceinterpolate (data is called functional there)
data2 = ft_checkdata(data, 'datatype', 'volume', 'inside', 'logical', 'feedback', 'yes', 'hasunits', 'yes');

% and we pretend we did some cluster statistics and have a label
data.negclusterslabelmat = nan(4,5,6);
data.negclusterslabelmat(1:10) = 1:10;

% this is how ft_checkdata is called from ft_sourceinterpolate (data is called functional there)
data3 = ft_checkdata(data, 'datatype', 'volume', 'inside', 'logical', 'feedback', 'yes', 'hasunits', 'yes');

