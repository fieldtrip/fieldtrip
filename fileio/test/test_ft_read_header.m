function test_ft_read_header

% TEST test_ft_read_header
% TEST ft_read_header

addpath('/home/common/matlab/fieldtrip/data/test/mfiles/');
datainfo = datasets;
datainfo = datainfo(1:17); % lfp does not yet work because of the data

for k = 1:numel(datainfo)
  subject  = datainfo(k);
  filename =  [subject.origdir,'original/',subject.type,subject.datatype,'/',subject.filename];

  % get header and event information
  if ~isempty(subject.dataformat)
    hdr   = ft_read_header(filename, 'headerformat', subject.dataformat);
  else
    hdr   = ft_read_header(filename);
  end
end
