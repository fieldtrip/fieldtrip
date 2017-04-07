function test_ft_preprocessing(datainfo, writeflag, version)

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_preprocessing ref_datasets

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

for k = 1:numel(datainfo)
  datanew = preprocessing10trials(datainfo(k), writeflag, version);

  fname = fullfile(datainfo(k).origdir,version,'raw',datainfo(k).type,['preproc_',datainfo(k).datatype]);
  load(fname);
  % these are per construction different if writeflag = 0;
  try, datanew = rmfield(datanew, 'cfg'); end
  try, data    = rmfield(data,    'cfg'); end
  % these can have subtle differences eg. in hdr.orig.FID
  data.hdr     = [];
  datanew2     = datanew; 
  datanew2.hdr = [];
  
  % do the comparison with the header removed, the output argument still contains the header
  %assert(isequaln(data, datanew2));  
  [ok,msg] = isalmostequal(data, datanew2,'reltol',eps*1e6);
  disp(['now you are in k=' num2str(k)]);
  if ~ok
    error('stored and computed data not identical: %s', msg{:});
  end
end


%----------------------------------------------------------
% subfunction to read in 10 trials of data
%----------------------------------------------------------
function [data] = preprocessing10trials(dataset, writeflag, version)

% --- HISTORICAL --- attempt forward compatibility with function handles
if ~exist('ft_preprocessing') && exist('preprocessing')
  eval('ft_preprocessing = @preprocessing;');
end
if ~exist('ft_read_header') && exist('read_header')
  eval('ft_read_header = @read_header;');
elseif ~exist('ft_read_header') && exist('read_fcdc_header')
  eval('ft_read_header = @read_fcdc_header;');
end
if ~exist('ft_read_event') && exist('read_event')
  eval('ft_read_event = @read_event;');
elseif ~exist('ft_read_event') && exist('read_fcdc_event')
  eval('ft_read_event = @read_fcdc_event;');
end

cfg            = [];
cfg.dataset    = fullfile(dataset.origdir,'original',dataset.type,dataset.datatype,'/',dataset.filename);
if writeflag,
  cfg.outputfile = fullfile(dataset.origdir,version,'raw',dataset.type,['preproc_',dataset.datatype '.mat']);
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

% JM 20170220: the neuromag306 raw.fif happens to have a bunch of zeros at
% the beginning of the file (up until ~sample 8300): this seems te be truly
% in the data, but creates problems if the generated data-structure is used
% downstream (e.g. in ft_sourceanalysis, with a Cf of all(zeros))), so for
% this type of data we shift the trl matrix a bit
if isfield(hdr, 'grad') && ft_senstype(hdr.grad, 'neuromag306')
  cfg.trl(:,1:2) = cfg.trl(:,1:2)+10000;
end

cfg.continuous = 'yes';
data           = ft_preprocessing(cfg);

if ~strcmp(version, 'latest') && str2num(version)<20100000
  % -- HISTORICAL --- older FieldTrip versions don't support outputfile
  save(cfg.outputfile, 'data');
end

