function data = avgoverdim(data, avgdim)

dimtok    = tokenize(data.dimord, '_');
avgdimnum = find(strcmp(avgdim, dimtok)); % the selected dimension as number
if length(avgdimnum)<1 && strcmp(avgdim, 'rpt'),
  avgdimnum = find(strcmp('rpttap', dimtok));
  avgdim    = 'rpttap';
end

if length(avgdimnum)<1 
  error('the "%s" dimension is not present in the data', avgdim)
elseif length(avgdimnum)>1
  error('cannot average over multiple dimensions at the same time')
end

reduceddim = dimlength(data);
if iscell(reduceddim)
  reduceddim = reduceddim{1}; % if reduceddim iscell, the next line of code cannot be executed
end
reduceddim(avgdimnum) = 1;

param = selparam(data);
for i=1:length(param)
  fprintf('averaging %s over %s\n', param{i}, avgdim);
  tmp = data.(param{i});
  tmp = reshape(nanmean(tmp, avgdimnum), reduceddim);
  data.(param{i}) = tmp;
end

switch avgdim
  case 'rpt'
    for i=1:length(param)
      fprintf('removing dimension %s from %s\n', avgdim, param{i});
      tmp = data.(param{i});
      tmp = reshape(tmp, [reduceddim(2:end) 1]);
      data.(param{i}) = tmp;
    end
    data.dimord = '';
    for i=2:length(dimtok)
      data.dimord = [data.dimord,'_',dimtok{i}];
    end
    data.dimord = data.dimord(2:end);
  case 'rpttap'
    for i=1:length(param)
      fprintf('removing dimension %s from %s\n', avgdim, param{i});
      warning('this is only allowed for cross-spectra and power-spectra');
      tmp = data.(param{i});
      tmp = reshape(tmp, [reduceddim(2:end) 1]);
      data.(param{i}) = tmp;
    end
    data.dimord = '';
    for i=2:length(dimtok)
      data.dimord = [data.dimord,'_',dimtok{i}];
    end
    data.dimord = data.dimord(2:end);
  
  case 'chan'
    data.label = avgoverlabel(data.label);
  case 'freq'
    data.freq = mean(data.freq);
  case 'time'
    data.time = mean(data.time);
  otherwise
    error('unknown dimension "%s"', avgdim);
end

if isfield(data, 'dim'),
  data.dim(avgdimnum) = [];
end
