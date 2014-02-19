function test_ft_read_sens

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_ft_read_sens
% TEST ft_read_sens

[dummy,ftpath] = ft_version();
addpath([ftpath '/test']);
datainfo = ref_datasets;

use      = match_str({datainfo.datatype},{'bti148' 'bti248' 'bti248grad' 'ctf151' 'ctf275' 'itab153' 'neuromag122' 'neuromag306'});
datainfo = datainfo(use);
for k = 1:numel(datainfo)
  subject  = datainfo(k);
  filename =  [subject.origdir,'original/',subject.type,'/',subject.datatype,'/',subject.filename];

  % get sensor information
  if ~isempty(subject.dataformat)
    sens = ft_read_sens(filename, 'headerformat', subject.dataformat);
  else
    sens = ft_read_sens(filename);
  end
end
