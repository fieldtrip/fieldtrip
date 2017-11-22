function test_ft_read_sens(datainfo, writeflag, version)

% MEM 1500mb
% WALLTIME 00:20:00

% TEST ft_read_sens

% writeflag determines whether the output should be saved to disk
% version determines the output directory

if nargin<1
  datainfo = ref_datasets;
end
if nargin<2
  writeflag = 0;
end
if nargin<3
  version = 'latest';
end

[dummy,ftpath] = ft_version();
addpath([ftpath '/test']);

% make a subselection of the MEG datasets
use      = match_str({datainfo.datatype},{'bti148' 'bti248' 'bti248grad' 'ctf151' 'ctf275' 'itab153' 'neuromag122' 'neuromag306'});
datainfo = datainfo(use);

for k = 1:numel(datainfo)
  dataset  = datainfo(k);
  filename =  [dataset.origdir,'original/',dataset.type,'/',dataset.datatype,'/',dataset.filename];
  
  % get sensor information
  if ~isempty(dataset.dataformat)
    sens = ft_read_sens(filename, 'headerformat', dataset.dataformat);
  else
    sens = ft_read_sens(filename);
  end
  
  if writeflag
    outputfile = fullfile(dataset.origdir,version,'sens',[dataset.senstype, '.mat']);
    save(outputfile, 'sens');
  end
  
end % for
