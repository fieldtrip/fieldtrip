function test_ft_checkdata


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


end


function sanityCheck(tmp)
  % sanity checks
  if ~isequal(size(tmp.sampleinfo), [4,2])
    error('sampleinfo is wrong');
  end

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
end

