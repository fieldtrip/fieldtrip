function data = seloverdim(data, seldim, sel, fb)

if nargin<4,
  fb = 1;
end

% get all XXXdimord fields
fn    = fieldnames(data);
selfn = find(~cellfun('isempty', strfind(fn, 'dimord')));
fn    = fn(selfn);
for k = 1:numel(fn)
  fndimord{k} = data.(fn{k});
end

% check which XXXdimord fields contain the seldim and keep only those
selx     = ~cellfun('isempty', strfind(fndimord, seldim));
fn       = fn(selx);
fndimord = fndimord(selx);

% extract the selected dimension as number
for k = 1:numel(fn)
  dimtok       = tokenize(fndimord{k}, '_');
  seldimnum{k} = find(strcmp(seldim, dimtok)); % the selected dimension as number
  if numel(seldimnum{k})<1 && strcmp(seldim, 'rpt'),
    %try rpttap
    seldimnum{k} = find(strcmp([seldim,'tap'], dimtok));
    tapflag      = 1;
  else
    tapflag      = 0;
  end
end

% if sum(~cellfun('isempty', seldimnum))<1
%   ft_error('the "%s" dimension is not present in the data', seldim)
% elseif any(cellfun(@numel, seldimnum)>1)
%   ft_error('cannot select over multiple dimensions at the same time')
% end

[reduceddim, fntmp] = dimlength(data);
selx       = find(ismember(fntmp, fn));
reduceddim = reduceddim(selx);
fntmp      = fntmp(selx);

% extract the fieldnames of the parameters of interest
if numel(fntmp)==1 && strcmp(fntmp{1}, 'dimord'),
  % data is not source data
  param      = selparam(data);
  reduceddim = repmat(reduceddim, [1 numel(param)]); 
  seldimnum  = repmat(seldimnum,  [1 numel(param)]);
else
  for k = 1:numel(fntmp)
    param{1,k} = fntmp{k}(1:end-6);
  end
end

% make the subselection
for i = 1:numel(param)
  if numel(seldimnum{i})==1,
    if fb, fprintf('selection %s along dimension %d\n', param{i}, seldimnum{i}); end
  else
    if fb, fprintf('selection %s along dimensions %d and %d\n', param{i}, seldimnum{i}(1), seldimnum{i}(2)); end
  end
    
  reduceddim{i}(seldimnum{i}) = numel(sel);
  tmp       = data.(param{i});
  iscelltmp = iscell(tmp);
  if ~iscelltmp,
    %temporarily convert to cell
    tmp = {tmp};
  else
    %keep cells but reduce seldimnum
    seldimnum{i}  = seldimnum{i} - 1; %FIXME this only works if cell is 1D
    reduceddim{i} = [reduceddim{i}(2:end) 1];
  end
 
  if seldimnum{i}>0,
    % subselection from each of the cells
    for j = 1:numel(tmp)
      if ~isempty(tmp{j}) && numel(seldimnum{i})==1,
        switch seldimnum{i}
          case 1
            tmp{j} = tmp{j}(sel,:,:,:,:,:,:,:,:);
          case 2
            tmp{j} = tmp{j}(:,sel,:,:,:,:,:,:,:);
          case 3
            tmp{j} = tmp{j}(:,:,sel,:,:,:,:,:,:);
          case 4
            tmp{j} = tmp{j}(:,:,:,sel,:,:,:,:,:);
          case 5
            tmp{j} = tmp{j}(:,:,:,:,sel,:,:,:,:);
          case 6
            tmp{j} = tmp{j}(:,:,:,:,:,sel,:,:,:);
          case 7
            tmp{j} = tmp{j}(:,:,:,:,:,:,sel,:,:);
          case 8
            tmp{j} = tmp{j}(:,:,:,:,:,:,:,sel,:);
          case 9
            tmp{j} = tmp{j}(:,:,:,:,:,:,:,:,sel);
          otherwise
            ft_error('the number of dimensions is too high');
        end
        tmp{j} = reshape(tmp{j}, reduceddim{i});
        
      elseif ~isempty(tmp{j}) && numel(seldimnum{i})>1,
        if all(seldimnum{i}==[1 2])
          tmp{j} = tmp{j}(sel,sel,:,:,:,:,:,:,:,:);
        elseif all(seldimnum{i}==[2 3])
          tmp{j} = tmp{j}(:,sel,sel,:,:,:,:,:,:,:);      
        else
          ft_error('selection of 2 dimensions simultaneously only works for [1 2] or [2 3] at present'); 
        end
        tmp{j} = reshape(tmp{j}, reduceddim{i});
      end
    end
  else
    % subselection of cells
    tmp = tmp(sel);
  end

  if ~iscelltmp,
    tmp = tmp{1};
  end
  data.(param{i}) = tmp;
end

switch seldim
  case 'rpt'
    if tapflag && isfield(data, 'cumtapcnt') && ~isempty(sel),
      sumtapcnt = cumsum([0;data.cumtapcnt(:)]);
      tapers    = zeros(1, sumtapcnt(end));
      for i=1:length(data.cumtapcnt)
        tapers(sumtapcnt(i)+1:sumtapcnt(i+1)) = i;
      end
      tmpsel    = [];
      tmpsel(1) = tapers(sel(1));
      for i=2:length(sel)
        if tapers(sel(i))~=tapers(sel(i-1))
	        tmpsel(end+1) = tapers(sel(i));
	      end
      end
    else
      tmpsel = sel;
    end
    if isfield(data, 'cumtapcnt'), data.cumtapcnt = data.cumtapcnt(tmpsel,:); end
    if isfield(data, 'cumsumcnt'), data.cumsumcnt = data.cumsumcnt(tmpsel,:); end
    
    % also try to adjust the trialinfo in the data
    if isfield(data, 'trialinfo')
      data.trialinfo = data.trialinfo(tmpsel, :);
    end
   
  case 'rpttap'
    % nothing to do
  case 'chan'
    data.label = data.label(sel);
  case 'freq'
    data.freq = data.freq(sel);
  case 'time'
    data.time = data.time(sel);
  case 'pos'
    data      = fixinside(data, 'logical');
    data.pos  = data.pos(sel,:);
    data.inside = data.inside(sel);
    data      = fixinside(data, 'index');
  otherwise
    ft_error('unknown dimension "%s"', seldim);
end
