function datanew = test_ft_preprocessing(datainfo, writeflag)

% TEST test_ft_preprocessing
% TEST ft_preprocessing test_datasets

% the optional writeflag determines whether the output 
% should be saved to disk

if nargin<2
  writeflag = 0;
end
if nargin<1
  datainfo = test_datasets;
end
for k = 1:numel(datainfo)
  datanew = preprocessing10trials(datainfo(k), writeflag);

  fname = [datainfo(k).origdir,'latest/raw/',datainfo(k).type,'preproc_',datainfo(k).datatype];
  load(fname);
  % these are per construction different if writeflag = 0;
  datanew = rmfield(datanew, 'cfg');
  data    = rmfield(data,    'cfg');
  % these can have subtle differences eg. in hdr.orig.FID
  data.hdr     = [];
  datanew2     = datanew; 
  datanew2.hdr = [];
  
  % do the comparison with the header removed, the output argument still
  % contains the header
  assert(isequalwithequalnans(data, datanew2));
end


%----------------------------------------------------------
% subfunction to read in 10 trials of data
%----------------------------------------------------------
function [data] = preprocessing10trials(dataset, writeflag)

cfg            = [];
cfg.dataset    = [dataset.origdir,'original/',dataset.type,dataset.datatype,'/',dataset.filename];
if writeflag,
  cfg.outputfile = [dataset.origdir,'latest/raw/',dataset.type,'preproc_',dataset.datatype];
end

% get header and event information
if ~isempty(dataset.dataformat)
  hdr   = ft_read_header(cfg.dataset, 'headerformat', dataset.dataformat);
  event = ft_read_event(cfg.dataset, 'eventformat', dataset.dataformat);
  
  cfg.dataformat   = dataset.dataformat;
  cfg.headerformat = dataset.dataformat;
else
  hdr   = ft_read_header(cfg.dataset);
  event = ft_read_event(cfg.dataset);
end

% create 10 1-second trials to be used as test-case 
begsample = ((1:10)-1)*round(hdr.Fs) + 1;
endsample = ((1:10)  )*round(hdr.Fs);
offset    = zeros(1,10);

cfg.trl   = [begsample(:) endsample(:) offset(:)];
sel       = cfg.trl(:,2)<=hdr.nSamples*hdr.nTrials;
cfg.trl   = cfg.trl(sel,:);

cfg.continuous = 'yes';
data           = ft_preprocessing(cfg);
