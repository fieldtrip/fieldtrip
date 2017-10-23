function data = leaveoneout(data)

dimtok = tokenize(data.dimord, '_');
rptdim = find(strcmp('rpt', dimtok)); % the selected dimension as number

if length(rptdim)<1
  ft_error('the ''rpt'' dimension is not present in the data');
elseif length(rptdim)>1
  ft_error('cannot jackknife over multiple dimensions at the same time');
elseif rptdim~=1
  ft_error('jackknife only works if replicates are in the first dimension of the data');
end

[reduceddim, fn] = dimlength(data);
if numel(fn)==1 && strcmp(fn{1}, 'dimord'),
  %data is not source data
  reduceddim = reduceddim{1};
else
  ft_error('jackknife not yet supported for source level data');
end
reduceddim(rptdim) = 1;

param = selparam(data);
for i=1:length(param)
  fprintf('computing jackknife %s\n', param{i});
  tmp    = data.(param{i});
  nrpt   = size(tmp, rptdim);
  sumtmp = reshape(nansum(tmp, rptdim), reduceddim);
  for k = 1:nrpt
    tmp(k,:,:,:,:,:,:) = (sumtmp - tmp(k,:,:,:,:,:,:))./(nrpt-1);
  end  
  data.(param{i}) = tmp;
end

data.method = 'jackknife';
