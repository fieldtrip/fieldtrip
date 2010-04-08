function data = seloverdim(data, seldim, sel)

dimtok    = tokenize(data.dimord, '_');
seldimnum = find(strcmp(seldim, dimtok)); % the selected dimension as number
if length(seldimnum)<1 && strcmp(seldim, 'rpt'),
  %try rpttap
  seldimnum = find(strcmp([seldim,'tap'], dimtok));
  tapflag   = 1;
else
  tapflag   = 0;
end

if length(seldimnum)<1
  error('the "%s" dimension is not present in the data', seldim)
elseif length(seldimnum)>1
  error('cannot select over multiple dimensions at the same time')
end

reduceddim = dimlength(data);
reduceddim(seldimnum) = length(sel);

param = selparam(data);
for i=1:length(param)
  fprintf('selection %s along dimension %d\n', param{i}, seldimnum);
  tmp = data.(param{i});
  switch seldimnum
    case 1
      tmp = tmp(sel,:,:,:,:,:,:,:,:);
    case 2
      tmp = tmp(:,sel,:,:,:,:,:,:,:);
    case 3
      tmp = tmp(:,:,sel,:,:,:,:,:,:);
    case 4
      tmp = tmp(:,:,:,sel,:,:,:,:,:);
    case 5
      tmp = tmp(:,:,:,:,sel,:,:,:,:);
    case 6
      tmp = tmp(:,:,:,:,:,sel,:,:,:);
    case 7
      tmp = tmp(:,:,:,:,:,:,sel,:,:);
    case 8
      tmp = tmp(:,:,:,:,:,:,:,sel,:);
    case 9
      tmp = tmp(:,:,:,:,:,:,:,:,sel);
    otherwise
      error('the number of dimensions is too high');
  end
  data.(param{i}) = reshape(tmp, reduceddim);
end

switch seldim
  case 'rpt'
    if tapflag && isfield(data, 'cumtapcnt'),,
      sumtapcnt = cumsum([0;data.cumtapcnt(:)]);
      tapers    = zeros(1, sumtapcnt(end));
      for i=1:length(data.cumtapcnt)
        tapers(sumtapcnt(i)+1:sumtapcnt(i+1)) = i;
      end
      tmpsel = unique(tapers(sel));
    else
      tmpsel = sel;
    end
    if isfield(data, 'cumtapcnt'), data.cumtapcnt = data.cumtapcnt(tmpsel); end
    if isfield(data, 'cumsumcnt'), data.cumsumcnt = data.cumsumcnt(tmpsel); end
   
    % also try to adjust the trl description in the configuration
    if isfield(data, 'cfg'), %try to locate the trl in the nested configuration
      trl = findcfg(data.cfg, 'trl');
    else
      trl = [];  
    end
    if isempty(trl) || size(trl,1)<length(tmpsel)
      % a trial definition is expected in each continuous data set
      warning('could not locate the correct trial definition ''trl'' in the data structure');
    else
      data.cfg.trlold = trl;
      data.cfg.trl    = trl(tmpsel,:);
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
    error('unknown dimension "%s"', seldim);
end
