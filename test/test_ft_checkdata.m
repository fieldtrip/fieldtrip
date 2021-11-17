function test_ft_checkdata

% MEM 2gb
% WALLTIME 00:20:00
% DEPENDENCY ft_checkdata

%% converting raw data to timelock data

% make some raw data with unequal time-axis, excluding 0
data = [];
data.label = {'1', '2'};

for m=[eps exp(1) pi 1:20]
  for n=.1:.1:m
    data.time{1} = -.5:(n/m):-.1;
    data.time{2} = -.5:(n/m):-.1;
    fsample = mean(diff(data.time{1}));
    if fsample <= 0 || isnan(fsample)
      continue;
    end
    for i=1:numel(data.time)
      data.trial{i} = rand(size(data.label, 2), size(data.time{i}, 2));
    end
    
    tmp = ft_checkdata(data, 'datatype', 'timelock', 'hassampleinfo', 'no');
    
    if (mean(diff(tmp.time)) - fsample > 12e-17)
      error('estimation of fsample does not match!')
    end
  end
end

for m=[eps exp(1) pi 1:20]
  for n=.1:.1:m
    data.time{1} = .1:(n/m):.5;
    data.time{2} = .1:(n/m):.5;
    fsample = mean(diff(data.time{1}));
    if fsample <= 0 || isnan(fsample)
      continue;
    end
    for i=1:numel(data.time)
      data.trial{i} = rand(size(data.label, 2), size(data.time{i}, 2));
    end
    
    tmp = ft_checkdata(data, 'datatype', 'timelock', 'hassampleinfo', 'no');
    
    if (mean(diff(tmp.time)) - fsample > 12e-17)
      error('estimation of fsample does not match!')
    end
  end
end

% make some raw data with strange time-axis
data = [];
data.label = {'1', '2'};

for m=[eps exp(1) pi 1:20]
  for n=.1:.1:m
    data.time{1} = [-(n.^2/m) -(n/m)];
    data.time{2} = [-(n.^2/m) -(n/m)];
    fsample = mean(diff(data.time{1}));
    if fsample <= 0 || isnan(fsample)
      continue;
    end
    for i=1:numel(data.time)
      data.trial{i} = rand(size(data.label, 2), size(data.time{i}, 2));
    end
    
    tmp = ft_checkdata(data, 'datatype', 'timelock', 'hassampleinfo', 'no');
    
    if (mean(diff(tmp.time)) - fsample > 12e-17)
      error('estimation of fsample does not match!')
    end
  end
  for n=eps^1.1:eps^1.1:eps
    data.time{1} = [-(n.^2/m) -(n/m)];
    data.time{2} = [-(n.^2/m) -(n/m)];
    fsample = mean(diff(data.time{1}));
    if fsample <= 0 || isnan(fsample)
      continue;
    end
    for i=1:numel(data.time)
      data.trial{i} = rand(size(data.label, 2), size(data.time{i}, 2));
    end
    
    tmp = ft_checkdata(data, 'datatype', 'timelock', 'hassampleinfo', 'no');
    
    if (mean(diff(tmp.time)) - fsample > 12e-17)
      error('estimation of fsample does not match!')
    end
  end
end

for m=[eps exp(1) pi 1:20]
  for n=eps:1:m
    data.time{1} = [-m -n];
    data.time{2} = [-m -n];
    fsample = mean(diff(data.time{1}));
    if fsample <= 0 || isnan(fsample)
      continue;
    end
    for i=1:numel(data.time)
      data.trial{i} = rand(size(data.label, 2), size(data.time{i}, 2));
    end
    
    if (n==eps)
      try
        tmp = ft_checkdata(data, 'datatype', 'timelock', 'hassampleinfo', 'no');
        if (mean(diff(tmp.time)) - fsample > 12e-17)
          warning('estimation of fsample does not match, but we''re near eps!')
        end
      catch
        warning('checkdata crashed, but we''re near eps!')
      end
    else
      tmp = ft_checkdata(data, 'datatype', 'timelock');
      if (mean(diff(tmp.time)) - fsample > 12e-17)
        error('estimation of fsample does not match!')
      end
    end
  end
  for n=eps^1.1:eps^1.1:eps
    data.time{1} = [-(n.^2/m) -(n/m)];
    data.time{2} = [-(n.^2/m) -(n/m)];
    fsample = mean(diff(data.time{1}));
    if fsample <= 0 || isnan(fsample)
      continue;
    end
    for i=1:numel(data.time)
      data.trial{i} = rand(size(data.label, 2), size(data.time{i}, 2));
    end
    
    tmp = ft_checkdata(data, 'datatype', 'timelock', 'hassampleinfo', 'no');
    
    if (mean(diff(tmp.time)) - fsample > 12e-17)
      error('estimation of fsample does not match!')
    end
  end
end

% make some raw data with unequal time-axis, excluding 0
data = [];
data.label = {'1', '2'};
data.time{1} = [-1.5 -1.28];
data.time{2} = [2.68 2.9];

for i=1:numel(data.time)
  data.trial{i} = rand(size(data.label, 2), size(data.time{i}, 2));
end

tmp = ft_checkdata(data, 'datatype', 'timelock', 'hassampleinfo', 'no');

if sum(tmp.time-[-1.5:0.22:3]) > 12e-17 || numel(tmp.time) ~= 21
  error('time axis is wrong');
end

% make some raw data with unequal time-axis, including 0 implicitly
data = [];
data.label = {'1', '2'};
data.time{1} = [-2 -1];
data.time{2} = [3 4];

for i=1:numel(data.time)
  data.trial{i} = rand(size(data.label, 2), size(data.time{i}, 2));
end

tmp = ft_checkdata(data, 'datatype', 'timelock', 'hassampleinfo', 'no');

if ~isequal(tmp.time, [-2:4])
  error('time axis is wrong');
end

% make some raw data with unequal time-axis, strictly < 0
% see bug 1477
data = [];
data.label = {'1', '2'};
data.time{1} = [-5 -4];
data.time{2} = [-3 -2];

for i=1:numel(data.time)
  data.trial{i} = rand(size(data.label, 2), size(data.time{i}, 2));
end

tmp = ft_checkdata(data, 'datatype', 'timelock', 'hassampleinfo', 'no');

if ~isequal(tmp.time, [-5:-2])
  error('time axis is wrong');
end

% make some raw data with unequal time-axis, strictly > 0
% related to bug 1477
data = [];
data.label = {'1', '2'};
data.time{1} = [4 5];
data.time{2} = [2 3];

for i=1:numel(data.time)
  data.trial{i} = rand(size(data.label, 2), size(data.time{i}, 2));
end

tmp = ft_checkdata(data, 'datatype', 'timelock', 'hassampleinfo', 'no');

if ~isequal(tmp.time, [2:5])
  error('time axis is wrong');
end

% make some raw data with unequal time-axis, including 0, with some jitter

data = [];
data.label = {'1', '2'};
data.time{1} = -.5:.25:1;
data.time{2} = -.25:.25:.25;
data.time{3} = .25:.25:1;
data.time{4} = -.5:.25:-.25;

for i=1:numel(data.time)
  data.time{i} = data.time{i} + (rand-0.5)/1000;
  data.trial{i} = rand(size(data.label, 2), size(data.time{i}, 2));
end
tmp = ft_checkdata(data, 'datatype', 'timelock', 'hassampleinfo', 'no');

%% converting raw data to timelock data

% make some raw data with unequal time-axis, including 0
data = [];
data.label = {'1', '2'};
data.time{1} = -.5:.25:1;
data.time{2} = -.25:.25:.25;
data.time{3} = .25:.25:1;
data.time{4} = -.5:.25:-.25;

for i=1:numel(data.time)
  data.trial{i} = rand(size(data.label, 2), size(data.time{i}, 2));
end

tmp = ft_checkdata(data, 'datatype', 'timelock');
sanityCheck(tmp);

%% shift time axis to be strictly positive
for i=1:numel(data.time)
  data.time{i} = data.time{i} + 0.6;
end

tmp = ft_checkdata(data, 'datatype', 'timelock');
sanityCheck(tmp);

%% shift time axis to be strictly negative
for i=1:numel(data.time)
  data.time{i} = data.time{i} - .6 - 1.1;
end

tmp = ft_checkdata(data, 'datatype', 'timelock');
sanityCheck(tmp);

%% make time-axis incredibly small
data.time{1} = -.5:.25:1;
data.time{2} = -.25:.25:.25;
data.time{3} = .25:.25:1;
data.time{4} = -.5:.25:-.25;
for i=1:numel(data.time)
  data.time{i} = (data.time{i}) ./ (10^-12);
end

tmp = ft_checkdata(data, 'datatype', 'timelock');
sanityCheck(tmp);

%% make time-axis awesomly huge
data.time{1} = -.5:.25:1;
data.time{2} = -.25:.25:.25;
data.time{3} = .25:.25:1;
data.time{4} = -.5:.25:-.25;
for i=1:numel(data.time)
  data.time{i} = data.time{i} .* (10^12);
end

tmp = ft_checkdata(data, 'datatype', 'timelock');
sanityCheck(tmp);

%% source to volume conversions, see https://github.com/fieldtrip/fieldtrip/pull/1920

% prepare the regular grid positions
dim = [4 5 6];
transform = eye(4);
[X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
pos = ft_warp_apply(transform, [X(:) Y(:) Z(:)]);

data = [];
data.dim = dim;
data.pos = pos;
data.time = 1:7;
data.freq = 1:8;

% this should be ignored in the conversion
data.mom = cell(prod(dim), 1);
for i=1:prod(dim)
  data.mom{i} = randn(3, length(data.time));
end
assert(strcmp(getdimord(data, 'mom'), '{pos}_ori_time'))

% this should be ignored in the conversion
data.filter = cell(prod(dim), 1);
data.filterdimord = '{pos}_chan_ori'; % otherwise it is detected as '{pos}_unknown_ori'
for i=1:prod(dim)
  data.filter{i} = randn(64, 3);
end
assert(strcmp(getdimord(data, 'filter'), '{pos}_chan_ori'))

data.coh = rand(prod(dim),1);
assert(strcmp(getdimord(data, 'coh'), 'pos'))
tmp = ft_checkdata(data, 'datatype', 'volume')
assert(strcmp(getdimord(tmp, 'coh'), 'dim1_dim2_dim3'))

data.coh = rand(prod(dim), length(data.time));
assert(strcmp(getdimord(data, 'coh'), 'pos_time'))
tmp = ft_checkdata(data, 'datatype', 'volume')
assert(strcmp(getdimord(tmp, 'coh'), 'dim1_dim2_dim3_time'))

data.coh = rand(prod(dim), length(data.time), length(data.freq));
assert(strcmp(getdimord(data, 'coh'), 'pos_time_freq'))
tmp = ft_checkdata(data, 'datatype', 'volume')
assert(strcmp(getdimord(tmp, 'coh'), 'dim1_dim2_dim3_time_freq'))

data.coh = rand(prod(dim), prod(dim));
assert(strcmp(getdimord(data, 'coh'), 'pos_pos'))
tmp = ft_checkdata(data, 'datatype', 'volume')
assert(strcmp(getdimord(tmp, 'coh'), 'dim1_dim2_dim3_dim1_dim2_dim3'))

data.coh = rand(prod(dim), prod(dim), length(data.time));
assert(strcmp(getdimord(data, 'coh'), 'pos_pos_time'))
tmp = ft_checkdata(data, 'datatype', 'volume')
assert(strcmp(getdimord(tmp, 'coh'), 'dim1_dim2_dim3_dim1_dim2_dim3_time'))

data.coh = rand(prod(dim), prod(dim), length(data.time), length(data.freq));
assert(strcmp(getdimord(data, 'coh'), 'pos_pos_time_freq'))
tmp = ft_checkdata(data, 'datatype', 'volume')
assert(strcmp(getdimord(tmp, 'coh'), 'dim1_dim2_dim3_dim1_dim2_dim3_time_freq'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sanityCheck(tmp)

% This does not need to be checked, because by construction
% a timelock structure should not contain a sampleinfo field
%if ~isequal(size(tmp.sampleinfo), [4,2])
%  error('sampleinfo is wrong');
%end

if ~isequal(tmp.time, [-.5:.25:1]) && ...
    ~isequal(tmp.time, [-.5:.25:1]+ 0.6) && ...
    ~isequal(tmp.time, [-.5:.25:1]- 1.1) &&  ...
    ~isequal(tmp.time, [-.5:.25:1]./ 10^-12) && ...
    ~isequal(tmp.time, [-.5:.25:1].* 10^-12)
  error('time axis is wrong');
end

% check individual trials
% note that we handle two channels here
if any(isnan(tmp.trial(1, :))) ...
    || any(isnan(tmp.trial(2, 3:8)))  || any(~isnan(tmp.trial(2, [1 2 9:14]))) ...
    || any(isnan(tmp.trial(3, 7:14))) || any(~isnan(tmp.trial(3, [1:6]))) ...
    || any(isnan(tmp.trial(4, 1:4)))  || any(~isnan(tmp.trial(4, [5:14])))
  error('nans are misplaced in .trial');
end


