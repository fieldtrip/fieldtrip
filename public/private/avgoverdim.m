function data = avgoverdim(data, avgdim)

% get all XXXdimord fields
fn    = fieldnames(data);
selfn = find(~cellfun('isempty', strfind(fn, 'dimord')));
fn    = fn(selfn);
for k = 1:numel(fn)
  fndimord{k} = data.(fn{k});
end

% check which XXXdimord fields contain the avgdim and keep only those
selx     = ~cellfun('isempty', strfind(fndimord, avgdim));
fn       = fn(selx);
fndimord = fndimord(selx);

% extract the selected dimension as number
for k = 1:numel(fn)
  dimtok       = tokenize(fndimord{k}, '_');
  avgdimnum{k} = find(strcmp(avgdim, dimtok)); % the selected dimension as number
  if numel(avgdimnum{k})<1 && strcmp(avgdim, 'rpt'),
    %try 'rpttap'
    avgdimnum{k} = find(strcmp('rpttap', dimtok));
    avgdim       = 'rpttap';
  end
end

if sum(~cellfun('isempty', avgdimnum))<1 
  error('the "%s" dimension is not present in the data', avgdim)
elseif any(cellfun(@numel, avgdimnum)>1)
  error('cannot average over multiple dimensions at the same time')
end

[reduceddim, fntmp] = dimlength(data);
selx       = find(ismember(fntmp, fn));
reduceddim = reduceddim(selx);
fntmp      = fntmp(selx);

% extract the fieldnames of the parameters of interest
if numel(fntmp)==1 && strcmp(fntmp{1}, 'dimord'),
  % data is not source data
  param      = selparam(data);
  reduceddim = repmat(reduceddim, [1 numel(param)]); 
  avgdimnum  = repmat(avgdimnum,  [1 numel(param)]);
else
  for k = 1:numel(fntmp)
    param{1,k} = fntmp{k}(1:end-6);
  end
end

for i = 1:numel(param)
  fprintf('averaging %s over %s\n', param{i}, avgdim);
  
  reduceddim{i}(avgdimnum{i}) = 1;
  tmp       = data.(param{i});
  iscelltmp = iscell(tmp);
  if ~iscelltmp,
    %temporarily convert to cell
    tmp = {tmp};
  else
    %keep cells but reduce avgdimnum
    avgdimnum{i}  = avgdimnum{i} - 1; %FIXME this only works if cell is 1D
    reduceddim{i} = [reduceddim{i}(2:end) 1];
  end 
  
  if avgdimnum{i}>0,
    % average each of the cells
    for j = 1:numel(tmp)
      tmp{j} = reshape(nanmean(tmp{j}, avgdimnum{i}), reduceddim{i}); 
    end
  else
    % not yet implemented
    error('averaging across cells is not yet possible');
  end
   
  if ~iscelltmp,
    tmp = tmp{1};
  end
  data.(param{i}) = tmp;
end

switch avgdim
  case 'rpt'
    for i = 1:length(param)
      fprintf('removing dimension %s from %s\n', avgdim, param{i});
      tmp = data.(param{i});
      tmp = reshape(tmp, [reduceddim{i}(2:end) 1]);
      data.(param{i}) = tmp;
    end
    data.dimord = '';
    for i = 2:length(dimtok)
      data.dimord = [data.dimord,'_',dimtok{i}];
    end
    data.dimord = data.dimord(2:end);
    if isfield(data, 'cumsumcnt'), data.cumsumcnt = sum(data.cumsumcnt); end
    if isfield(data, 'cumtapcnt'), data.cumtapcnt = sum(data.cumtapcnt); end
  
  case 'rpttap'
    for i=1:length(param)
      fprintf('removing dimension %s from %s\n', avgdim, param{i});
      warning('this is only allowed for cross-spectra and power-spectra');
      tmp = data.(param{i});
      tmp = reshape(tmp, [reduceddim{i}(2:end) 1]);
      data.(param{i}) = tmp;
    end
    data.dimord = '';
    for i=2:length(dimtok)
      data.dimord = [data.dimord,'_',dimtok{i}];
    end
    data.dimord = data.dimord(2:end);
    if isfield(data, 'cumsumcnt'), data.cumsumcnt = sum(data.cumsumcnt); end
    if isfield(data, 'cumtapcnt'), data.cumtapcnt = sum(data.cumtapcnt); end
  
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
  data.dim(avgdimnum{1}) = [];
end
