function test_ft_read_header

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_read_header

[ftver, ftpath] = ft_version;
addpath([ftpath '/test']);
datainfo = ref_datasets;
datainfo = datainfo(1:17); % lfp does not yet work because of the data

for k = 1:numel(datainfo)
  subject  = datainfo(k);
  filename =  [subject.origdir,'original/',subject.type,'/',subject.datatype,'/',subject.filename];

  % get header and event information
  if ~isempty(subject.dataformat)
    hdr   = ft_read_header(filename, 'headerformat', subject.dataformat);
  else
    hdr   = ft_read_header(filename);
  end
end
